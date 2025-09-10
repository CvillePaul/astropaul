from glob import glob
from pathlib import Path
import re
from sqlite3 import Connection

from astropaul.database import database_path
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import numpy as np
import pandas as pd

ptab_data_line = r"""(?x)
    (?P<ra>\d{5})
    (?P<dec>[\+-]\d{4})\s+
    (?P<filename>\S+)\s+
    (?P<byear>\d\d.\d\d\d\d)\s+
    (?P<red_chisq>[\d\.]+)\s+
    (?P<seeing>[\d\.]+)\s+
    (?P<position_angle>[\d\.]+)\s+
    (?P<separation>[\d\.]+)\s+
    (?P<delta_mag>[\d\.]+)\s+
    (?P<filter>\d+)\s*
    (?P<notes>.*)
"""


def process_dssi_ptab_file(
    filename: str, target_names: list[str], target_coords: SkyCoord, observations: pd.DataFrame
) -> pd.DataFrame:
    match_column = "start_jd" # compare this column to the byear in the ptab file
    use_column = "mid_jd" # when match is found, use this column for the output file
    observation_times = observations.copy()[["target_name", match_column, use_column, "mid_utc", "wavelengths"]]
    observation_times.reset_index(inplace=True)
    match_threshold = TimeDelta(61 * u.min) # byear with 4 decimals should be 52 min, but needed 61 to get all lines to match
    answer = pd.DataFrame(
        columns=[
            "Target Name",
            "Mid JD",
            "Mid UTC",
            "Reduced ChiSq",
            "Seeing",
            "Position Angle",
            "Separation",
            "Delta Mag",
            "Filter",
            "Filename",
            "Notes",
        ]
    )

    with open(filename) as f:
        for line in f.readlines():
            if not (match := re.match(ptab_data_line, line)):
                print(f"Skipped: {line}", end="")
                continue
            fields = match.groupdict()
            # use coordinates in file & coords of known targets as crosscheck of stated target id
            ra = f"{fields["ra"][:2]}:{fields["ra"][2:4]}:{int(fields["ra"][4:5])/60:.0f}"
            dec = f"{fields["dec"][:3]}:{fields["dec"][3:5]}:00"
            result_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            idx, sep2d, _ = result_coord.match_to_catalog_sky(target_coords, nthneighbor=1)
            if sep2d[0] > 0.5 * u.deg:
                print(f"Skipped unmatched target: {result_coord}")
                continue  # skip unknown targets
            nearest_known_target_name = target_names[idx]
            if "target_id" in fields and not fields["target_id"] in nearest_known_target_name: # comparing this way avoids assuming all targets are TICs
                raise ValueError(f"Target ID {fields["target_id"]} has coordinates closest to {nearest_known_target_name}")
            # replace the less precise byear with the more precise jd from the matching observation row
            byear = fields["byear"]
            result_time = Time(2000 + float(byear), format="byear")
            target_criteria = observation_times["target_name"] == nearest_known_target_name
            time_criteria = np.abs((Time(observation_times[match_column], format="jd") - result_time).jd) < match_threshold.jd
            filter = fields["filter"]
            filter_criteria = np.array([filter in wavelengths for wavelengths in observation_times["wavelengths"]])
            matching_observations = observation_times[target_criteria & time_criteria & filter_criteria]
            if matching_observations.empty:
                raise ValueError(f"No observation found for target {nearest_known_target_name} and byear {byear}")
            observation_time = Time(matching_observations.iloc[0][use_column], format="jd")
            remaining_wavelengths = matching_observations.iloc[0]["wavelengths"].replace(filter, "")
            observation_times.loc[matching_observations.iloc[0]["index"], "wavelengths"] = remaining_wavelengths
            # write out results
            answer.loc[len(answer)] = [
                nearest_known_target_name,
                observation_time.jd,
                observation_time.iso[:19],
                fields["red_chisq"],
                fields["seeing"],
                fields["position_angle"],
                fields["separation"],
                fields["delta_mag"],
                fields["filter"],
                fields["filename"],
                fields["notes"],
            ]
        return answer


def process_dssi_ptab_files(files: list[str], out_dir: str = ".", database: str = None, verbose: bool = False) -> None:
    if not database:
        database = database_path()
    with Connection(database) as conn:
        targets = pd.read_sql("select * from targets;", conn)
        target_names = targets["target_name"]
        target_coords = SkyCoord(targets["ra"], targets["dec"], unit="deg")
        observations = pd.read_sql("select * from dssi_observations;", conn)
    out_path = Path(out_dir)
    for file_pattern in files:
        for specific_file in glob(file_pattern):
            print(f"Processing {specific_file}")
            dssi_results = process_dssi_ptab_file(specific_file, target_names, target_coords, observations)
            dssi_results = dssi_results.sort_values(["Target Name", "Mid UTC"])
            out_file = out_path / Path(specific_file).name.replace(".ptab", ".csv")
            dssi_results.to_csv(out_file, index=False)
            if verbose:
                print(f"Wrote {len(dssi_results)} result lines to {out_file}")
