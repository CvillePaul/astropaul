import argparse
from collections import Counter
from datetime import datetime, timedelta
from glob import glob
import os
import re
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

wavelengths_optical = "692, 880"
wavelengths_ir = "1450"


def get_file_contents(file: str) -> tuple[datetime, list[str]]:
    with open(file) as f:
        first_line = f.readline()
        date_string = first_line.replace("DSSI OBSERVING at the ARC Telescope, ", "")
        observing_date = datetime.strptime(date_string, "%d %b %Y UT\n")
        contents = f.readlines()
    return observing_date, contents


def find_observation_lines(lines: list[str]) -> tuple[list[str], list[str]]:
    """Take list of lines, split into tuple of observation lines and non-observation lines."""
    observation_lines = []
    non_observation_lines = []
    for line in lines:
        if line.count(":") >= 4:  # all observation lines follow this pattern
            observation_lines.append(line)
        else:
            non_observation_lines.append(line)
    return observation_lines, non_observation_lines


observation_line_patterns = [
    (
        "Standard Pattern",
        r"""(?x)
            (?P<target_name>.{7,13})\s+
            (?P<image_num>\d{1,3})\s+
            (?P<hours>\d\d):(?P<minutes>\d\d)\s+
            (?P<gain_1>\d{1,3})\s+
            (?P<gain_2>\d{1,3})\s+
            (?P<ra>\d\d:\d\d:\d\d\.\d+)\s+
            (?P<dec>[+|-]{0,1}\d\d:\d\d:\d\d\.\d+)\s+
            (?P<pmra>[0-9.+-]*)\s+
            (?P<pmdec>[0-9.+-]*)\s+
            (?P<mag>[0-9.-]+)\s*
            (?P<notes>.*)""",
    ),
    # (
    #     "No Image Num, Gains, or Time",
    #     r"""(?x)
    #         (?P<target_name>\"{0,1}.{7,13}\"{0,1})\s+
    #         (?P<ra>\d\d:\d\d:\d\d\.\d+)\s+
    #         (?P<dec>[+|-]{0,1}\d\d:\d\d:\d\d\.\d+)\s+
    #         (?P<pmra>[0-9.+-]*)\s+
    #         (?P<pmdec>[0-9.+-]*)\s+
    #         (?P<mag>[0-9.-]+)\s*
    #         (?P<notes>.*)""",
    # ),
    (
        "Infrared Observations",
        r"""(?x)
            (?P<target_name>.{7})\s+
            (?P<image_beg>\d{1,3})-(?P<image_end>\d{1,3})\s+
            (?P<image_ir>\d{1,3})\s+
            (?P<hours>\d\d):(?P<minutes>\d\d)\s+
            (?P<gain>\d{1,3})\s+
            (?P<ra>\d\d:\d\d:\d\d\.\d+)\s+
            (?P<dec>[+|-]{0,1}\d\d:\d\d:\d\d\.\d+)\s+
            (?P<pmra>[0-9\.+-]*)\s+
            (?P<pmdec>[0-9\.+-]*)\s+
            (?P<mag>[0-9\.]+)\s*
            (?P<notes>.*)""",
    ),
    (
        "Single Gain Value",
        r"""(?x)
            (?P<target_name>.{7,15})\s+
            (?P<image_num>\d{1,4})\s+
            (?P<gain>\d{1,3})\s+
            (?P<hours>\d\d):(?P<minutes>\d\d)\s+
            (?P<ra>\d\d:\d\d:\d\d\.\d+)\s+
            (?P<dec>[+|-]{0,1}\d\d:\d\d:\d\d\.\d+)\s+
            (?P<pmra>[0-9\.+-]*)\s+
            (?P<pmdec>[0-9\.+-]*)\s+
            (?P<mag>[0-9\.]+)\s*
            (?P<notes>.*)""",
    ),
]


def process_observation_lines(observation_lines: list[str], utc_date: datetime) -> tuple[Table, Counter, list[str]]:
    """Extract data from observation line, return Table of observations, count of types of line, and failed lines."""

    dssi_sequences = Table(
        names=[
            "Target Name",
            "Wavelengths",
            "Image Number",
            "UTC DateTime",
            "Time JD",
            "Gain 1",
            "Gain 2",
            "RA",
            "Dec",
            "PMRA",
            "PMDec",
            "Mag",
            "Notes",
        ],
        dtype=[
            "str",
            "str",
            "int",
            "str",
            "float",
            "int",
            "int",
            "float",
            "float",
            "float",
            "float",
            "float",
            "str",
        ],
    )

    observation_types = Counter()
    failed_lines = []

    for line_num, line in enumerate(observation_lines):
        fields = {}
        for pattern_name, regex_pattern in observation_line_patterns:
            if match := re.match(regex_pattern, line):
                observation_types.update(pattern_name)
                break
        if match:
            fields = {**match.groupdict(), **fields}
            try:
                if not "pmra" in fields or fields["pmra"] == "":
                    fields["pmra"] = 0
                    fields["pmdec"] = 0
            except:
                pass
            if "gain" in fields:  # assume gains are the same for both arms if only one gain specified
                fields["gain_1"] = fields["gain"]
                fields["gain_2"] = fields["gain"]
            if "hours" in fields:
                obs_time = Time(utc_date + timedelta(hours=int(fields["hours"]), minutes=int(fields["minutes"])))
                datetime_utc = str(obs_time.utc)
                datetime_jd = obs_time.jd
            else:
                datetime_utc = str(utc_date)
                datetime_jd = 0
            if "image_ir" in fields:
                observations = [
                    (image_num, wavelengths_optical)
                    for image_num in range(int(fields["image_beg"]), int(fields["image_end"]) + 1)
                ]
                observations.append((fields["image_ir"], wavelengths_ir))
            else:
                observations = [(fields["image_num"], wavelengths_optical)]
            coord = SkyCoord(ra=fields["ra"], dec=fields["dec"], unit=(u.hourangle, u.deg))
            for image_num, wavelengths in observations:
                try:
                    dssi_sequences.add_row(
                        [
                            fields["target_name"].replace('"', ""),
                            wavelengths,
                            image_num,
                            datetime_utc,
                            datetime_jd,
                            fields["gain_1"],
                            fields["gain_2"],
                            coord.ra.value,
                            coord.dec.value,
                            float(fields["pmra"]),
                            float(fields["pmdec"]),
                            fields["mag"],
                            fields["notes"],
                        ]
                    )
                except Exception as e:
                    failed_lines.append(f"{line_num:4d}: {line.strip()}")
        else:
            failed_lines.append(f"{line_num:4d}: {line.strip()}")
    return dssi_sequences, observation_types, failed_lines


def determine_dssi_observations(dssi_sequences: Table) -> Table:
    # add a column to the sequences table that indicates observation number
    prev_target = ""
    speckle_session = 0
    dssi_sequences["DSSI Session"] = 0
    dssi_sequences.sort("Time JD")
    for sequence in dssi_sequences:
        target_name = sequence["Target Name"]
        if target_name != prev_target:
            speckle_session += 1
            prev_target = target_name
        sequence["DSSI Session"] = speckle_session
    # group all sequences for a given observation, and summarize them in a new table
    dssi_observations = Table(
        names=[
            "Target Name",
            "DSSI Session",
            "Start JD",
            "Mid JD",
            "End JD",
            "Mid UTC",
            "Num Sequences",
            "Wavelengths",
        ],
        dtype=[str, str, float, float, float, str, int, str],
    )

    obs_by_session = dssi_sequences.group_by(["DSSI Session", "Target Name", "Wavelengths"])
    for keys, sequences in zip(obs_by_session.groups.keys, obs_by_session.groups):
        start_time = sequences["Time JD"].min()
        end_time = sequences["Time JD"].max()
        mid_time = (end_time + start_time) / 2
        mid_utc = str(Time(mid_time, format="jd").iso)[:19] if mid_time > 0 else ""
        dssi_observations.add_row(
            (
                keys["Target Name"],
                str(keys["DSSI Session"]),
                start_time,
                mid_time,
                end_time,
                mid_utc,
                len(sequences),
                keys["Wavelengths"],
            )
        )

    return dssi_observations


def parse_olist_file(file: str) -> tuple[datetime, Table, Counter, list[str], list[str]]:
    observing_date, contents = get_file_contents(file)
    observation_lines, non_observation_lines = find_observation_lines(contents)
    dssi_sequences, observation_types, failed_lines = process_observation_lines(observation_lines, observing_date)
    if dssi_sequences and len(dssi_sequences) > 0:
        dssi_observations = determine_dssi_observations(dssi_sequences)
    else:
        dssi_observations = None
    return observing_date, dssi_observations, observation_types, failed_lines, non_observation_lines


def parse_olist_files(files: list[str], out_dir: str = ".", verbose: bool = False) -> None:
    overall_observation_types = Counter()
    overall_failed_lines = {}
    for file_pattern in files:
        for specific_file in glob(file_pattern):
            observation_date, dssi_observations, observation_types, failed_lines, _ = parse_olist_file(specific_file)
            if not dssi_observations or len(dssi_observations) == 0:
                if verbose:
                    print(f"{0:4d} sequences, {len(failed_lines):4d} failed lines from {specific_file}")
                continue
            sequence_file = f"DSSI Observations {observation_date:%Y-%m-%d}.csv"
            sequence_file = os.path.join(out_dir, sequence_file)
            dssi_observations.write(sequence_file, overwrite=True)
            overall_observation_types += observation_types
            overall_failed_lines[specific_file] = failed_lines
            if verbose:
                print(f"{len(dssi_observations):4d} sequences, {len(failed_lines):4d} failed lines from {specific_file}")
    if verbose:
        print("Processing complete")
        print(f"{sum(overall_observation_types.values())} sequences found")
        for observation_type, count in overall_observation_types.items():
            print(f"    {count:4d}: {observation_type}")
        if (num_failed := sum(map(len, overall_failed_lines.values()))) > 0:
            print(f"{num_failed} observation lines failed matches:")
            for file, lines in overall_failed_lines.items():
                if len(lines) > 0:
                    print(f"{file}")
                    for line in lines:
                        print(f"    {line}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="parse_olist_files",
        description="Find lines in olist file(s) representing dssi sequences, collect them into csv file of same name",
    )
    parser.add_argument("olist_files", nargs="+", help="List of .olist files to process")
    parser.add_argument("-o", "--out_dir", default=".", help="Directory for resulting csv files (default = %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output stats about file processing")
    args = parser.parse_args()
    parse_olist_files(args.olist_files, out_dir=args.out_dir, verbose=args.verbose)
