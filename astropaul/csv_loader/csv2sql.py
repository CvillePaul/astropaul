import argparse
import os

import pandas as pd
import sqlalchemy as sa 

def convert_naming_style(string:str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return string.lower().replace(" ", "_")

def infer_sqlalchemy_type(series:pd.Series) -> str:
    if pd.api.types.is_integer_dtype(series):
        return sa.Integer
    elif pd.api.types.is_float_dtype(series):
        return sa.Float
    else:
        return sa.String  # Default to string for categorical/text data

def validate_csv_schemas(csv_files) -> dict[str, object]:
    reference_df = pd.read_csv(csv_files[0])  # Use first file as reference
    reference_columns = list(reference_df.columns)
    reference_types = {col: infer_sqlalchemy_type(reference_df[col]) for col in reference_columns}

    for csv_file in csv_files[1:]:
        df = pd.read_csv(csv_file)
        if list(df.columns) != reference_columns:
            raise ValueError(f"Schema mismatch in {csv_file}. Expected columns: {reference_columns}, Found: {list(df.columns)}")
        for col in df.columns:
            if infer_sqlalchemy_type(df[col]) != reference_types[col]:
                raise ValueError(f"Data type mismatch in {csv_file} for column {col}.")
    return reference_types

def create_table(engine, metadata, table_name:str, column_map:dict[str, object]):
    table_columns = []
    # table_columns.append(sa.Column("id", sa.Integer, primary_key=True, autoincrement=True))
    for column_name, column_type in column_map.items():
        table_columns.append(sa.Column(convert_naming_style(column_name), column_type))
    table = sa.Table(table_name, metadata, *table_columns)
    return table

def find_foreign_keys(engine, metadata):
    """Identify and apply foreign key relationships based on column name equality."""
    table_columns = {table.name: set(col.name for col in table.columns) for table in metadata.tables.values()}

    for table_name, columns in table_columns.items():
        table = metadata.tables[table_name]
        for col in columns:
            for other_table, other_columns in table_columns.items():
                if table_name != other_table and col in other_columns:
                    fk = sa.ForeignKey(f"{other_table}.{col}")
                    fk_col_name = f"fk_{other_table}"
                    table.append_column(sa.Column(fk_col_name, sa.Integer, fk))
                    print(f"Foreign key added: {table_name}.{fk_col_name} -> {other_table}.{col}")

def insert_csv_data(engine:sa.Engine, csv_files:list[str], table_name:str) -> int:
    num_rows = 0
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df.columns = [convert_naming_style(column) for column in df.columns]
        df.to_sql(table_name, engine, if_exists="append", index=False)
        num_rows += len(df)
    return num_rows

def csv2sql(base_dir:str, outfile:str, verbose:bool = False):
    try:
        if os.path.exists(outfile):
            os.remove(outfile)
        engine = sa.create_engine(f"sqlite:///{outfile}")
        metadata = sa.MetaData()
        csv_files_by_table = {}
        for table_dir in os.listdir(base_dir):
            table_path = os.path.join(base_dir, table_dir)
            if not os.path.isdir(table_path):
                continue  # Skip files, only process directories
            table_name = convert_naming_style(table_dir)
            csv_files = [
                os.path.join(table_path, f) for f in os.listdir(table_path) if f.lower().endswith(".csv") and f[0] != "#"
            ]
            csv_files_by_table[table_name] = csv_files
            if not csv_files:
                if verbose:
                    print(f"Skipping empty directory {table_dir}")
                continue  # Skip empty directories
            column_map = validate_csv_schemas(csv_files)
            create_table(engine, metadata, table_name, column_map)
        metadata.create_all(engine)
        # load data into the new tables
        for table_name, csv_files in csv_files_by_table.items():
            num_rows = insert_csv_data(engine, csv_files, table_name)
            if verbose:
                print(f"{num_rows:5d} rows inserted into table {table_name}")
    finally:
        engine.dispose()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="csv2sql",
        description="",
    )
    parser.add_argument("-o", "--outfile", default="astropaul.db", help="Name of SQLite3 file to create (default = %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output file info as each file is processed")
    parser.add_argument("directory", default=".", help="Base directory of CSV files (default = %(default)s")
    args = parser.parse_args()
    csv2sql(args.directory, outfile=args.outfile, verbose=args.verbose)
