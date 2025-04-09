import argparse
import os
import re

import networkx as nx
import pandas as pd
import sqlalchemy as sa


def convert_naming_style(string: str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return string.lower().replace(" ", "_")


def infer_sqlalchemy_type(series: pd.Series) -> str:
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
    reference_types = {convert_naming_style(column): column_type for column, column_type in reference_types.items()}
    return reference_types


def create_table(engine, metadata, table_name: str, column_map: dict[str, object]):
    table_columns = []
    # table_columns.append(sa.Column("id", sa.Integer, primary_key=True, autoincrement=True))
    for column_name, column_type in column_map.items():
        table_columns.append(sa.Column(column_name, column_type))
    table = sa.Table(table_name, metadata, *table_columns)
    return table


def find_foreign_keys(tables: set[str]) -> nx.DiGraph:
    """Identify and apply foreign key relationships based on column name equality."""
    foreign_key_re = r"(?P<child_table>.*)\.(?P<child_column>.*) -> (?P<source_table>.*)\.(?P<source_column>.*)"
    dependencies = {}  # child_table: {child_column: (source_table, source_column)}
    with open("foreign_keys.txt") as f:
        for line_num, line in enumerate(f.readlines()):
            if line.startswith("#"):
                continue  # skip comment lines
            if not (match := re.match(foreign_key_re, line)):
                raise ValueError(f"Syntax error on line {line_num}: {line}")
            child_table = convert_naming_style(match.groupdict()["child_table"])
            child_column = convert_naming_style(match.groupdict()["child_column"])
            child_table_dependencies = dependencies.get(child_table, {})
            if child_column in child_table_dependencies:
                raise ValueError(f"Child table.column {child_table}.{child_column} appear multiple times in foreign keys list")
            child_table_dependencies[child_column] = (
                convert_naming_style(match.groupdict()["source_table"]),
                convert_naming_style(match.groupdict()["source_column"]),
            )
            dependencies[child_table] = child_table_dependencies
    edges = [
        (source_table, child_table)
        for child_table, child_table_dependencies in dependencies.items()
        for source_table, _ in child_table_dependencies.values()
    ]
    dependency_graph = nx.DiGraph()
    for table in tables:
        dependency_graph.add_node(table)
    for child_table, child_table_deps in dependencies.items():
        for _, (source_table, _) in child_table_deps.items():
            dependency_graph.add_edge(source_table, child_table)
    if not nx.is_directed_acyclic_graph(dependency_graph):
        raise ValueError(f"Dependency graph is not acyclic.  Dependencies are: {dependencies}")
    table_order = list(nx.topological_sort(dependency_graph))
    return table_order, dependencies


def insert_csv_data(
    engine: sa.Engine, metadata: sa.MetaData, csv_files: list[str], table_name: str, fk_constraints: dict[str, tuple[str, str]]
) -> int:
    validation_values = {}
    for child_column, (source_table_name, source_column) in fk_constraints.items():
        source_table = metadata.tables[source_table_name]
        with engine.connect() as conn:
            result = conn.execute(sa.select(source_table.c[source_column]))
            validation_values[child_column] = set([row[0] for row in result])
    num_rows = 0
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df.columns = [convert_naming_style(column) for column in df.columns]
        for child_column in fk_constraints.keys():
            for value in df[child_column]:
                if not value in validation_values[child_column]:
                    raise ValueError(f"Value {value} not found in table {source_table} column {source_column}")
        df.to_sql(table_name, engine, if_exists="append", index=False)
        num_rows += len(df)
    return num_rows


def csv2sql(base_dir: str, outfile: str, verbose: bool = False):
    engine = None
    try:
        if os.path.exists(outfile):
            os.remove(outfile)
        engine = sa.create_engine(f"sqlite:///{outfile}")
        metadata = sa.MetaData()
        csv_files_by_table = {}
        table_dirs = {}  # table_name: table_dir
        for table_dir in os.listdir(base_dir):
            table_path = os.path.join(base_dir, table_dir)
            if not os.path.isdir(table_path):
                continue  # Skip files, only process directories
            table_name = convert_naming_style(table_dir)
            table_dirs[table_name] = table_dir
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
        table_order, fk_constraints = find_foreign_keys(set(table_dirs.keys()))
        # load data into the new tables
        for table_name in table_order:
            table_name_sql = convert_naming_style(table_name)
            csv_files = csv_files_by_table[table_name_sql]
            num_rows = insert_csv_data(engine, metadata, csv_files, table_name_sql, fk_constraints.get(table_name, {}))
            if verbose:
                print(f"{num_rows:5d} rows inserted into table {table_name}")
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
