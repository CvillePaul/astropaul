from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import pandas as pd

from astropaul.database import database_connection

filter_specs = {
    "y": 517,
    "I": 824,
}

def process_tokovinin_email_catalog(file: Path, data_dir: Path) -> None:
    # read in known targets
    with database_connection() as conn:
        targets = pd.read_sql("select target_name, ra, dec from targets;", conn)
    targets_coord = SkyCoord(targets["ra"], targets["dec"], unit=u.deg)

    # tk = pd.read_table(file, usecols=[], index_col=False)
    tk = pd.read_fwf(
        file, 
        colspecs=[
            (0, 10),
            (11, 26),
            (28, 42),
            (44, 51),
            (53, 54),
            (56, 57),
            (58, 59),
            (62, 68),
            (70, 75),
            (79, 82),
            (83, 84),
            (85, 93),
        ],
    )
    tk_coord = SkyCoord(
        [f"{coord[0:2]}h{coord[2:4]}.{coord[4:5]}m{coord[5:8]}d{coord[8:10]}m" for coord in tk["WDS"]], unit=(u.hourangle, u.deg)
    )
    tk_time = Time([2000 + dyear for dyear in tk["YR-2000"]], format="decimalyear")
    idx, dist, _ = tk_coord.match_to_catalog_sky(targets_coord)

    # use some columns to make table of observations
    tk_obs = pd.DataFrame()
    tk_obs["Target Name"] = targets.iloc[idx]["target_name"]
    tk_obs["UTC"] = tk_time.iso
    tk_obs["JD"] = tk_time.jd
    for col in ["F", "N", "B", "X"]:
        tk_obs[col] = [val for val in tk[col]]
    tk_obs.to_csv(data_dir / "HRCam Observations" / file.name.replace(".txt", ".csv"), index=False)

    # use some columns to make table of speckle detections
    tk_det = pd.DataFrame()
    tk_det["Target Name"] = targets.iloc[idx]["target_name"]
    tk_det["Mid UTC"] = tk_time.iso
    tk_det["Observation ID"] = [f"{utc}~{target}" for target, utc in zip(tk_det["Target Name"], tk_det["Mid UTC"])]
    tk_det["Instrument"] = "HRCam"
    tk_det["Mid JD"] = tk_time.jd
    for old_col, new_col in [("Rho", "Separation"), ("Theta", "Position Angle"), ("dm", "Delta Mag")]:
        tk_det[new_col] = [val for val in tk[old_col]]
    tk_det["Reduced ChiSq"] = -1.0 # float("nan")
    tk_det["Seeing"] = -1.0 # float("nan")
    tk_det["Filter"] = [filter_specs[filter] for filter in tk["F"]]
    tk_det["Filename"] = " "
    tk_det["Notes"] = " "
    tk_det["B"] = [val for val in tk["B"]]
    tk_det = tk_det[tk_det["B"] == "B"].drop("B", axis=1)
    tk_det.to_csv(data_dir / "Speckle Detections" / file.name.replace(".txt", ".csv"), index=False, na_rep="NaN")
