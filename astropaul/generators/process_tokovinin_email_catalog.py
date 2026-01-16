from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import pandas as pd

from astropaul.database import database_connection

def process_tokovinin_email_catalog():
    # read in known targets
    with database_connection() as conn:
        targets = pd.read_sql("select target_name, ra, dec from targets;", conn)
    targets_coord = SkyCoord(targets["ra"], targets["dec"], unit=u.deg)

    tk = pd.read_table("../../Files/HRCam/From Tokovenin 2025-04-14.txt", sep=r"\s+", index_col=False)

    tk_coord = SkyCoord(
        [f"{coord[0:2]}h{coord[2:4]}.{coord[4:5]}m{coord[5:8]}d{coord[8:10]}m" for coord in tk["WDS"]], unit=(u.hourangle, u.deg)
    )
    tk_time = Time([2000 + dyear for dyear in tk["YR-2000"]], format="decimalyear")
    idx, dist, _ = tk_coord.match_to_catalog_sky(targets_coord)
    targets.iloc[idx]

    # use some columns to make table of observations
    tk_obs = pd.DataFrame()
    tk_obs["Target Name"] = targets.iloc[idx]["target_name"]
    tk_obs["UTC"] = tk_time.iso
    tk_obs["JD"] = tk_time.jd
    for col in ["F", "N", "B"]:
        tk_obs[col] = [val for val in tk[col]]
    tk_obs.to_csv("../../Data/HRCam Observations/From Tokovenin 2025-04-14.csv", index=False)

    # use some columns to make table of speckle detections
    tk_det = pd.DataFrame()
    tk_det["Target Name"] = targets.iloc[idx]["target_name"]
    tk_det["Mid UTC"] = tk_time.iso
    tk_det["Mid JD"] = tk_time.jd
    for old_col, new_col in [("Rho", "Separation"), ("Theta", "Position Angle"), ("dm", "Delta Mag")]:
        tk_det[new_col] = [val for val in tk[old_col]]
    tk_det["Reduced ChiSq"] = float("nan")
    tk_det["Seeing"] = float("nan")
    tk_det["Filter"] = 0
    tk_det["Filename"] = " "
    tk_det["Notes"] = " "
    tk_det.to_csv("../../Data/Speckle Detections/From Tokovenin 2025-04-14.csv", index=False)
