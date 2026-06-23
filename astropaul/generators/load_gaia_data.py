import argparse
import os
from pathlib import Path
import sys

from astroquery.gaia import Gaia
from astropy.table import Column
import pandas as pd

from astropaul.database import database_connection, database_path

def load_gaia_data(out_dir: Path = None, verbose: bool = True) -> None:
    """Make a csv of Gaia DR3 data for all possible entries in the TESS database table."""

    out_dir = out_dir or database_path().parent

    # get the gaia DR2 source IDs for all targets that have TESS data
    dr2_query = """
        select cast(gaia as INT) gaia_dr2, target_name
        from tess_ticv8
        where gaia is not null
        ;
    """
    with database_connection() as conn:
        dr2_to_target_name = pd.read_sql(dr2_query, conn, index_col="gaia_dr2").to_dict()["target_name"]

    # retrieve dr3 data using the dr2 id list
    gaia_query = f"""
        select *
        FROM gaiaedr3.gaia_source AS gaia
        JOIN gaiaedr3.dr2_neighbourhood AS xmatch
          ON gaia.source_id = xmatch.dr3_source_id
        WHERE xmatch.dr2_source_id IN ({",".join(map(str, dr2_to_target_name.keys()))})
        """
    try:
        job = Gaia.launch_job(gaia_query)
        results = job.get_results()
    except Exception as e:
        print(e)

    results.add_column(Column([dr2_to_target_name[dr2] for dr2 in results["dr2_source_id"]], name="Target Name"), index=0)
    out_file = out_dir / "Gaia DR3" / "Gaia DR3.csv"
    results.write(out_file, overwrite=True)

    out_file = out_dir / "Catalog Members" / "Gaia DR3.csv"
    members = pd.DataFrame()
    members["Target Name"] = dr2_to_target_name.values()
    members["Catalog Name"] = "Gaia DR3"
    members["Catalog ID"] = dr2_to_target_name.keys()
    members.to_csv(out_file, index=False)

    if verbose:
        print(f"Wrote {len(results)} rows to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description="Download Gaia DR3 catalog data for every target with a TESS TICv8 catalog entry.",
    )
    parser.add_argument("-d", "--database", required=True, help="Database file to use as input")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory to put resulting CSV file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output stats about file processing")
    args = parser.parse_args()
    load_gaia_data(args.database, args.out_dir, verbose=args.verbose)
