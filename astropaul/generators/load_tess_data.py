import argparse
import os
from pathlib import Path
import sys

from astropaul.database import database_connection
from astroquery.mast import Catalogs
import pandas as pd


def load_tess_data(database: str = None, out_dir: str = ".", verbose: bool = False) -> None:
    # retrieve all objects by TIC ID from the TIC catalog
    # TIC fields: https://mast.stsci.edu/api/v0/_t_i_cfields.html
    catalog_query = "select target_name, catalog_id from catalog_members where catalog_name = 'TESS TICv8';"
    with database_connection(database) as conn:
        target_ids = pd.read_sql(catalog_query, conn)
    tic_table = Catalogs.query_criteria(catalog="Tic", ID=target_ids["catalog_id"]).to_pandas()
    unique_ids = set(target_ids["catalog_id"])
    if len(tic_table) != len(unique_ids):
        raise ValueError(f"Tic query returned {len(tic_table)} instead of {len(unique_ids)} entries")
    tess_data = target_ids.merge(tic_table, left_on="catalog_id", right_on="ID")
    tess_data.rename(columns={"target_name": "Target Name", "ID": "TIC", "ra": "TESS RA", "dec": "TESS Dec"}, inplace=True)
    tess_data.drop(columns=["catalog_id"], inplace=True)
    # now reorder columns for prettiness
    first_columns = ["Target Name", "TIC", "Vmag", "Teff"]
    second_columns = [col for col in tess_data.columns if not col in first_columns]
    tess_data = tess_data[first_columns + second_columns]
    out_file = Path(out_dir) / "TESS TICv8.csv"
    tess_data.to_csv(out_file, index=False)
    if verbose:
        print(f"Wrote {len(tess_data)} rows to {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description="Download TESS catalog data for every target with a TESS TICv8 catalog entry.",
    )
    parser.add_argument("-d", "--database", required=True, help="Database file to use as input")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory to put resulting CSV file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output stats about file processing")
    args = parser.parse_args()
    load_tess_data(args.database, args.out_dir, verbose=args.verbose)
