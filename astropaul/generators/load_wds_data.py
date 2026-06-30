from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Column, QTable
import astropy.units as u
import pandas as pd

from astropaul.database import base_path, database_connection, database_path


def load_wds_data(out_dir: Path = None, wds_file: Path = None, verbose: bool = True) -> None:

    out_dir = out_dir or database_path().parent
    wds_file = wds_file or base_path() / "Files" / "Washington Double Star Catalog" / "wdsweb_summ2.ecsv"

    target_query = """
        select target_name, ra, dec 
        from targets
        where target_type in ("QuadEB", "SextEB")
        ;"""
    with database_connection() as conn:
        targets = pd.read_sql(target_query, conn)
    target_coords = SkyCoord(ra=targets["ra"], dec=targets["dec"], unit=(u.deg, u.deg))

    wds = QTable.read(wds_file)
    wds_coords = SkyCoord(ra=wds["RA"], dec=wds["Dec"], unit=(u.deg, u.deg))

    # find matches for each target in the WDS
    target_names = [""] * len(wds)
    idx, d2d, _ = target_coords.match_to_catalog_sky(wds_coords)
    for i, (wds_index, distance) in enumerate(zip(idx, d2d)):
        if distance < 1 * u.arcmin:
            target_names[wds_index] = targets.loc[i, "target_name"]
    wds.rename_column("ID", "WDS ID")
    wds.add_column(Column(target_names, name="Target Name"), index=0)
    wds = wds[wds["Target Name"] != ""]

    # write out WDS ID mapping
    out_file = out_dir / "Catalog Members" / "WDS.csv"
    catalog_members = pd.DataFrame()
    catalog_members["Target Name"] = wds["Target Name"]
    catalog_members["Catalog Name"] = "Washington Double Star"
    catalog_members["Catalog ID"] = wds["WDS ID"]
    catalog_members.to_csv(out_file, index=False)

    # write out matched WDS entries
    out_file = out_dir / "WDS" / "WDS.csv"
    wds.write(out_file, overwrite=True)


    if verbose:
        print(f"{len(targets)} targets had {len(wds)} WDS matches, written to {out_file}")
