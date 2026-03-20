import argparse
import os
from pathlib import Path
import sys

from astropaul.database import database_connection
from astropy.table import Table
from astroquery.utils.tap.core import TapPlus
import pandas as pd


def create_gaia_catalog_members(database:str = None, out_dir: str = ".", verbose: bool = False) -> None:
    """Create Gaia DR3 associations from target names"""

def load_gaia_data(database: str = None, out_dir: str = ".", verbose: bool = False) -> None:
    """Make a csv of Gaia DR3 data for all possible entries in the TESS database table."""


    # get table of all our targets that are mapped to a TESS TIC id
    tic_query = "select target_name, catalog_id from catalog_members where catalog_name = 'TESS TICv8';"
    with database_connection(database) as conn:
        tic_targets = Table.from_pandas(pd.read_sql(tic_query, conn))

    tic_targets = tic_targets[:1]

    with TapPlus(url="https://mast.stsci.edu/api/v0/tap/tap") as mast_tap: 

        upload_name = "target_list"
        mast_tap.upload_table(upload_resource=upload_name, table=tic_targets)

        # Build ADQL query
        mast_query = f"""
        SELECT targets.target_name, dr3.*
        FROM tap_upload.{upload_name} targets
        JOIN tic.tic as tess
        ON tess.ticid = targets.catalog_id
        JOIN gaia.dr3lite as dr3
        ON tess.dr3_source_id = dr3.source_id
        """

        job = mast_tap.launch_job(mast_query)
        results = job.get_results()

        mast_tap.delete_table(upload_name)

    out_file = Path(out_dir) / "Gaia DR3.csv"
    results.write(out_file, overwrite=True)
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
