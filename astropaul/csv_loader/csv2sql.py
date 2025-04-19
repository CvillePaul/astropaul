import argparse
import configparser
import os
from pathlib import Path
import re

import networkx as nx
import pandas as pd
import sqlalchemy as sa


def db_naming_style(name: str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return name.lower().replace(" ", "_")


def series_to_column_type(series: pd.Series) -> str:
    if pd.api.types.is_integer_dtype(series):
        return sa.Integer
    elif pd.api.types.is_float_dtype(series):
        return sa.Float
    else:
        return sa.String  # Default to string for categorical/text data


def string_to_column_type(column_type_str: str) -> object:
    match column_type_str.lower():
        case "str":
            column_type = sa.String
        case "int":
            column_type = sa.Integer
        case "float":
            column_type = sa.Float
        case _:
            column_type = None
    return column_type


def validate_csv_schemas(
    data_files: list[tuple[Path, pd.DataFrame]], specified_columns: dict[str, object]
) -> dict[str, object]:
    reference_file, reference_df = data_files[0]  # Use first file as reference
    all_columns = set(reference_df.columns)
    column_types = {}
    for column, column_type in specified_columns.items():
        if column not in all_columns:
            raise ValueError(f"Column {column} specified in config file but not present in {reference_file}")
        column_types[column] = column_type
    inferred_columns = {
        column: series_to_column_type(reference_df[column])
        for column in reference_df.columns
        if column not in specified_columns
    }
    column_types = {**column_types, **inferred_columns}
    for file_name, data_file in data_files[1:]:  # validate all remaining files against schema we just built
        df_columns = set(data_file.columns)
        if df_columns != all_columns:
            raise ValueError(f"Schema mismatch in {file_name}. Expected columns: {all_columns}, Found: {df_columns}")
        for column, column_type in inferred_columns.items():
            this_type = series_to_column_type(data_file[column])
            if this_type != column_type:
                raise ValueError(f"Data type mismatch in {file_name} for column {column}: {this_type} instead of {column_type}")
    return column_types


def create_table(engine, metadata, table_name: str, column_map: dict[str, object]):
    table_columns = []
    # table_columns.append(sa.Column("id", sa.Integer, primary_key=True, autoincrement=True))
    for column_name, column_type in column_map.items():
        table_columns.append(sa.Column(db_naming_style(column_name), column_type))
    table = sa.Table(db_naming_style(table_name), metadata, *table_columns)
    return table


def determine_table_order(table_names: list[str], constraints: dict[str, dict[str, tuple[str, str]]]) -> list[str]:
    dependency_graph = nx.DiGraph()
    for table_name in table_names:
        dependency_graph.add_node(table_name)
    for child_table, child_table_deps in constraints.items():
        for _, (source_table, _) in child_table_deps.items():
            dependency_graph.add_edge(source_table, child_table)
    if not nx.is_directed_acyclic_graph(dependency_graph):
        raise ValueError(f"Dependency graph is not acyclic.  Dependencies are: {constraints}")
    table_order = list(nx.topological_sort(dependency_graph))
    return table_order


def insert_csv_data(
    engine: sa.Engine,
    metadata: sa.MetaData,
    data_files: list[tuple[str, pd.DataFrame]],
    table_name: str,
    constraints: dict[str, tuple[str, str]],
) -> int:
    validation_values = {}
    for child_column, (source_table_name, source_column) in constraints.items():
        source_table = metadata.tables[db_naming_style(source_table_name)]
        with engine.connect() as conn:
            result = conn.execute(sa.select(source_table.c[db_naming_style(source_column)]))
            validation_values[child_column] = set([row[0] for row in result])
    num_rows = 0
    for file_name, data_file in data_files:
        for child_column in constraints.keys():
            valid_column_values = validation_values[child_column]
            for value in data_file[child_column].unique():
                if not value in valid_column_values:
                    raise ValueError(
                        f"In file {file_name} for table {table_name}, value {value} not found in {source_table}.{source_column}"
                    )
        db_frame = data_file.copy()
        db_frame.columns = [db_naming_style(column) for column in db_frame.columns]
        db_frame.to_sql(db_naming_style(table_name), engine, if_exists="append", index=False)
        num_rows += len(data_file)
    return num_rows


def parse_table_config(table_dir: Path) -> tuple[dict[str, object], dict[str, str]]:
    config_file = table_dir / f"{table_dir.name}.ini"
    if config_file.exists():
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(config_file)
        specified_columns = {}
        if "columns" in config.sections():
            for column, column_type_str in config.items("columns"):
                if not (column_type := string_to_column_type(column_type_str)):
                    raise ValueError(f"Unknown column type {column_type_str} in file {config_file}")
                specified_columns[column] = column_type

        constraints = {}  # key is column name, value is (table name, column name) of foreign table
        if "constraints" in config.sections():
            for column, foreign_column in config.items("constraints"):
                parts = foreign_column.split(".")
                if not parts or len(parts) != 2:
                    raise ValueError(f"Bad constraint {foreign_column} for column {column} in {config_file}")
                foreign_table, foreign_column = parts
                constraints[column] = (foreign_table, foreign_column)
    else:
        specified_columns, constraints = {}, {}
    return specified_columns, constraints


def get_data_files(data_dir: Path) -> list[tuple[str, pd.DataFrame]]:
    data_files = []
    for file in data_dir.iterdir():
        if not file.name.lower().endswith("csv"):
            continue  # only process CSV files
        if file.name[0] == "#":
            continue  # skip files with "commented out" names
        data_files.append((file, pd.read_csv(file)))
    return data_files


def csv2sql(base_dir: str, outfile: str, verbose: bool = False):
    engine = None
    try:
        outfile_path = Path(outfile)
        if outfile_path.exists():
            outfile_path.unlink()
        engine = sa.create_engine(f"sqlite:///{outfile_path}")
        metadata = sa.MetaData()

        # find files by subdir, validate their schemas, and create empty tables
        data_files_by_table = {}  # key is table name, value is list of (filename, DataFrame)
        constraints = {}  # key is table name, value is dict of column: (table, column) of foreign table
        for table_dir in Path(base_dir).iterdir():
            table_name = table_dir.name
            if not table_dir.is_dir():
                continue  # Skip files, only process directories
            data_files = get_data_files(table_dir)
            if not data_files:
                if verbose:
                    print(f"Skipping directory with no data files: {table_dir}")
                continue
            data_files_by_table[table_name] = data_files
            specified_columns, table_constraints = parse_table_config(table_dir)
            constraints[table_name] = table_constraints
            column_map = validate_csv_schemas(data_files, specified_columns)
            create_table(engine, metadata, table_name, column_map)
            if verbose:
                print(f"Created table {table_name} with {len(column_map)} columns")
        metadata.create_all(engine)

        # load data into the new tables
        table_order = determine_table_order(data_files_by_table.keys(), constraints)
        for table_name in table_order:
            data_files = data_files_by_table[table_name]
            num_rows = insert_csv_data(engine, metadata, data_files, table_name, constraints.get(table_name, {}))
            if verbose:
                print(f"{num_rows:5d} rows inserted into table {table_name} from {len(data_files)} files")

        if verbose:
            print(f"Database {outfile_path} created")
    finally:
        if engine:
            engine.dispose()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="csv2sql",
        description="",
    )
    parser.add_argument(
        "-o", "--outfile", default="astropaul.db", help="Name of SQLite3 file to create (default = %(default)s)"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Output file info as each file is processed")
    parser.add_argument("directory", default=".", help="Base directory of CSV files (default = %(default)s")
    args = parser.parse_args()
    csv2sql(args.directory, outfile=args.outfile, verbose=args.verbose)
