import argparse
from glob import glob
import os
from pathlib import Path
from sqlite3 import Connection
import sys

from astroquery.mast import Catalogs
import pandas as pd


def load_tess_data(database: str, out_dir: str, verbose: bool = False) -> None:
    # retrieve all objects by TIC ID from the TIC catalog
    # TIC fields: https://mast.stsci.edu/api/v0/_t_i_cfields.html
    with Connection(database) as conn:
        id_list = pd.read_sql("select catalog_id from catalog_members where catalog_name = 'TESS TICv8';", conn)["catalog_id"]
    id_list = [id.replace("TIC", "").strip() for id in id_list]
    tic_table = Catalogs.query_criteria(catalog="Tic", ID=id_list)
    # tic_table.rename_column("ID", "Identifier")
    if len(tic_table) != len(id_list):
        raise ValueError(f"Tic query returned {len(tic_table)} instead of {len(id_list)} entries")
    out_file = Path(out_dir) / "TESS TICv8.csv"
    tic_table.write(out_file, overwrite=True)
    if verbose:
        print(f"Wrote {len(tic_table)} rows to {out_file}")


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
