from pathlib import Path

from astroquery.gaia import Gaia
from astropy.table import Column
import astropy.units as u
import pandas as pd

from astropaul.database import database_connection, database_path

def load_gaia_bailer_jones(out_dir: Path = None, verbose: bool = True) -> None:
    """For every DR3 id in Catalog Members, download Bailer-Jones distance information (https://iopscience.iop.org/article/10.3847/1538-3881/abd806)"""

    out_dir = out_dir or database_path().parent

    # get the gaia DR3 source IDs for all targets
    dr3_query = """
        select catalog_id, target_name
        from catalog_members
        where catalog_name = 'Gaia DR3'
        ;
    """
    with database_connection() as conn:
        dr3_to_target_name = pd.read_sql(dr3_query, conn, index_col="catalog_id").to_dict()["target_name"]

    # retrieve bailer-jones distance from the dr3 ids
    gaia_query = f"""
        select source_id, r_med_geo, r_hi_geo, r_lo_geo
        FROM external.gaiaedr3_distance
        WHERE source_id IN ({",".join(map(str, dr3_to_target_name.keys()))})
        """
    try:
        job = Gaia.launch_job(gaia_query)
        results = job.get_results()
    except Exception as e:
        print(e)

    # add target name as first column, remove source_id, rename distance column
    results.add_column(Column([dr3_to_target_name[str(dr3)] for dr3 in results["source_id"]], name="Target Name"), index=0)
    results.remove_column("source_id")
    results.rename_column("r_med_geo", "Bailer Jones Distance")
    results.rename_column("r_lo_geo", "Bailer Jones Lower Limit")
    results.rename_column("r_hi_geo", "Bailer Jones Upper Limit")

    # write out distance data
    out_file = out_dir / "Gaia Bailer Jones" / "Gaia Bailer Jones.csv"
    results.write(out_file, overwrite=True)

    if verbose:
        print(f"Wrote {len(results)} rows to {out_file}")

