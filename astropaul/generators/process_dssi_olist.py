import argparse
from collections import Counter
from datetime import datetime, timedelta
from glob import glob
import os
import re
from sqlite3 import connect

from astropaul.database import database_path
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import pandas as pd

wavelengths_optical = "692, 880"
wavelengths_ir = "1500"


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
        if line.count(":") > 4:  # all observation lines follow this pattern
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
    #     "No Image Num, Gains, or Time", # Jimmy said these lines shouldn't be used
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

    dssi_sequences = pd.DataFrame()
    observation_types = Counter()
    failed_lines = []

    for line_num, line in enumerate(observation_lines):
        fields = {}
        for pattern_name, regex_pattern in observation_line_patterns:
            if match := re.match(regex_pattern, line):
                observation_types.update([pattern_name])
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
            coord = SkyCoord(fields["ra"], fields["dec"], unit=(u.hourangle, u.deg))
            for image_num, wavelengths in observations:
                try:
                    dssi_sequences = pd.concat([
                        dssi_sequences, 
                        pd.DataFrame({
                            "Target Name": [fields["target_name"].replace('"', "")],
                            "Wavelengths": [wavelengths],
                            "Image Number": [image_num],
                            "UTC DateTime": [datetime_utc],
                            "Time JD": [datetime_jd],
                            "Gain 1": [fields["gain_1"]],
                            "Gain 2": [fields["gain_2"]],
                            "RA": [coord.ra.value],
                            "Dec": [coord.dec.value],
                            "PMRA": [float(fields["pmra"])],
                            "PMDec": [float(fields["pmdec"])],
                            "Mag": [fields["mag"]],
                            "Notes": [fields["notes"]],
                        })])
                except Exception as e:
                    failed_lines.append(f"{line_num:4d}: {line.strip()}")
        else:
            failed_lines.append(f"{line_num:4d}: {line.strip()}")
    return dssi_sequences, observation_types, failed_lines


def determine_dssi_observations(dssi_sequences: Table, starting_session_num: int) -> Table:
    # add a column to the sequences table that indicates observation number
    prev_target = ""
    speckle_session = starting_session_num
    session_numbers = []
    for index, row in dssi_sequences.sort_values("Time JD").iterrows():
        target_name = row["Target Name"]
        if target_name != prev_target:
            speckle_session += 1
            prev_target = target_name
        session_numbers.append(speckle_session)
    dssi_sequences["DSSI Session"] = session_numbers
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

    for keys, sequences in dssi_sequences.groupby(["DSSI Session", "Target Name", "Wavelengths"]):
        start_time = sequences["Time JD"].min()
        end_time = sequences["Time JD"].max()
        mid_time = (end_time + start_time) / 2
        mid_utc = str(Time(mid_time, format="jd").iso)[:19] if mid_time > 0 else ""
        dssi_observations.add_row(
            (
                keys[1],
                str(keys[0]),
                start_time,
                mid_time,
                end_time,
                mid_utc,
                len(sequences),
                keys[2],
            )
        )
    return dssi_observations.to_pandas().sort_values("Mid JD")


def parse_olist_file(file: str, known_targets: pd.DataFrame, starting_session_num: int) -> tuple[datetime, Table, Counter, list[str], list[str]]:
    observing_date, contents = get_file_contents(file)
    observation_lines, non_observation_lines = find_observation_lines(contents)
    dssi_sequences, observation_types, failed_lines = process_observation_lines(observation_lines, observing_date)
    if not dssi_sequences.empty:
        if not known_targets.empty:
            target_coords = SkyCoord(known_targets["ra"], known_targets["dec"], unit="deg")
            sequence_coords = SkyCoord(dssi_sequences["RA"], dssi_sequences["Dec"], unit="deg")
            idx, sep2d, _ = sequence_coords.match_to_catalog_sky(target_coords)
            matches = sep2d < 20 * u.arcsec
            dssi_sequences.loc[matches, "Target Name"] = known_targets.loc[idx]["target_name"][matches].to_list()
        dssi_observations = determine_dssi_observations(dssi_sequences, starting_session_num)
    else:
        dssi_observations = pd.DataFrame()
    return observing_date, dssi_observations, observation_types, failed_lines, non_observation_lines


def parse_olist_files(files: list[str], out_dir: str = ".", database: str = None, verbose: bool = False) -> None:
    overall_observation_types = Counter()
    overall_failed_lines = {}
    session_num = 0
    if not database:
        database = database_path()
    with connect(database) as conn:
        known_targets = pd.read_sql("select target_name, ra, dec from targets;", conn)
    for file_pattern in files:
        for specific_file in glob(file_pattern):
            observation_date, dssi_observations, observation_types, failed_lines, _ = parse_olist_file(specific_file, known_targets, session_num)
            session_num += len(dssi_observations)
            if dssi_observations.empty:
                if verbose:
                    print(f"{0:4d} sequences, {len(failed_lines):4d} failed lines from {specific_file}")
                continue
            sequence_file = f"DSSI Observations {observation_date:%Y-%m-%d}.csv"
            sequence_file = os.path.join(out_dir, sequence_file)
            dssi_observations.to_csv(sequence_file, index=False)
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
    parser.add_argument("-d", "--database", default=None, help="Database with known targets for use in name substitution")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output stats about file processing")
    args = parser.parse_args()
    parse_olist_files(args.olist_files, out_dir=args.out_dir, database=args.database, verbose=args.verbose)
