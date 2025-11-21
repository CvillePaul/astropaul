import argparse
from pathlib import Path

import networkx as nx
import pandas as pd
import sqlalchemy as sa

from .db_utils import *
from .table_config import *
from .data_transformations import *

def create_table(metadata: sa.MetaData, table_config: TableConfig, table_metadata: pd.DataFrame) -> int:
    table_columns = []
    primary_key_type = table_config.key_generator.get_key_type()
    table_columns.append(sa.Column("id", primary_key_type, primary_key=True))
    for transformation in table_config.transformations:
        for sql_name, sql_type in transformation.get_sql_columns():
            table_columns.append(sa.Column(sql_name, sql_type))
    sa.Table(table_config.table_name, metadata, *table_columns)
    # write info from ini file about units to the metadata table
    for column_name, unit in table_config.units.items():
        table_metadata.loc[len(table_metadata)] = [table_config.table_name, column_name, "unit", unit]
    # write metadata from the ini file to the metadata table
    if table_config.metadata:
        for metadata_field, metadata_value in table_config.metadata.items():
            table_metadata.loc[len(table_metadata)] = [table_config.table_name, None, metadata_field, metadata_value]
    return len(table_columns)


def determine_table_order(table_configs: dict[str, TableConfig]) -> list[str]:
    dependency_graph = nx.DiGraph()
    for table_name in table_configs.keys():
        dependency_graph.add_node(table_name)
    for child_table, table_config in table_configs.items():
        for _, (source_table, _) in table_config.constraints.items():
            dependency_graph.add_edge(source_table, child_table)
    if not nx.is_directed_acyclic_graph(dependency_graph):
        raise ValueError(f"Dependency graph is not acyclic.  Dependencies are: {table_config.constraints}")
    table_order = list(nx.topological_sort(dependency_graph))
    return table_order


def insert_csv_data(
    engine: sa.Engine,
    metadata: sa.MetaData,
    data_files: list[tuple[str, pd.DataFrame]],
    table_config: TableConfig,
) -> int:
    # retrieve values of foreign key constraints for validation before inserting data
    validation_values = {}
    for child_column, (source_table_name, source_column) in table_config.constraints.items():
        source_table = metadata.tables[source_table_name]
        with engine.connect() as conn:
            result = conn.execute(sa.select(source_table.c[source_column]))
            validation_values[child_column] = set([row[0] for row in result])
    # validate and insert rows from each data file
    num_rows = 0
    num_constraint_violations = None
    for file_name, df_from_csv in data_files:
        # verify constrained columns have allowed values
        for child_column in table_config.constraints.keys():
            constraint_violations = []
            valid_column_values = validation_values[child_column]
            for value in df_from_csv[child_column].unique():
                if not value in valid_column_values:
                    constraint_violations.append(value)
            if len(constraint_violations) == 0:
                continue
            if not num_constraint_violations:
                num_constraint_violations = 0
            num_constraint_violations += len(constraint_violations)
            message = f"In file {file_name}, column {child_column} contains invalid values {constraint_violations}"
            if "log" in table_config.constraint_options:
                print(message)
            if "skip" in table_config.constraint_options:
                violation_mask = df_from_csv[child_column].isin(constraint_violations)
                df_from_csv = df_from_csv[~violation_mask]
            else:
                raise ValueError(message)
        if df_from_csv.empty:
            continue
        # now build the DataFrame we'll convert to sql
        df_for_sql = pd.DataFrame()
        df_for_sql["id"] = table_config.key_generator.generate_primary_keys(df_from_csv)
        for transformation in table_config.transformations:
            df_for_sql = transformation.do_transformation(df_from_csv, df_for_sql)
        df_for_sql.to_sql(table_config.table_name, engine, if_exists="append", index=False)
        num_rows += len(df_from_csv)
    return num_rows, num_constraint_violations


def get_data_files(data_dir: Path) -> list[tuple[str, pd.DataFrame]]:
    data_files = []
    for file in data_dir.iterdir():
        if not file.name.lower().endswith("csv"):
            continue  # only process CSV files
        if file.name[0] == "#":
            continue  # skip files with "commented out" names
        df = pd.read_csv(file)
        df.columns = [string_to_db_style(column_name) for column_name in df.columns]
        data_files.append((file, df))
    return data_files


def csv2sql(base_dir: str, outfile: str=None, verbose: bool=False):
    engine = None
    try:
        outfile_path = Path(outfile) if outfile else database_path()
        if outfile_path.exists():
            outfile_path.unlink()
        engine = sa.create_engine(f"sqlite:///{outfile_path}")
        metadata = sa.MetaData()
        table_metadata = pd.DataFrame(columns=["table_name", "column_name", "value_type", "value"])

        # find files by subdir, validate their schemas, and create empty tables
        data_files_by_table = {}  # key is table name, value is list of (filename, DataFrame)
        table_configs = {}  # key is table name, value is TableConfig
        for table_dir in Path(base_dir).iterdir():
            table_name = table_dir.name
            if not table_dir.is_dir():
                continue  # Skip files, only process directories
            data_files = get_data_files(table_dir)
            if len(data_files) == 0:
                if verbose:
                    print(f"Skipping directory with no data files found: {table_dir}")
                continue
            data_files_by_table[string_to_db_style(table_name)] = data_files
            config_file = table_dir / f"{table_name}.ini"
            table_config = TableConfig(config_file, data_files[0][1])
            table_configs[table_config.table_name] = table_config
            num_columns = create_table(metadata, table_config, table_metadata)
            if verbose:
                print(f"Created table {table_name} with {num_columns} columns, {len(table_config.constraints)} constraints")
        metadata.create_all(engine)
        table_metadata.to_sql("table_metadata", engine, index=False)
        if verbose:
            print(f"Created metadata table with {len(table_metadata)} rows")

        # load data into the new tables
        table_order = determine_table_order(table_configs)
        for table_name in table_order:
            table_config = table_configs[table_name]
            if verbose:
                print(f"Table {db_style_to_string(table_name)}: ", end="")
            data_files = data_files_by_table[table_name]
            num_rows, num_violations = insert_csv_data(engine, metadata, data_files, table_config)
            if verbose:
                if num_violations:
                    voilation_text = f", {num_violations} rows skipped"
                else:
                    voilation_text = ""
                print(f"{num_rows} rows inserted from {len(data_files)} files{voilation_text}")

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
    parser.add_argument("-d", "--directory", default=".", help="Base directory of CSV files (default = %(default)s)")
    parser.add_argument(
        "-o", "--outfile", default="astropaul.db", help="Name of SQLite3 file to create (default = %(default)s)"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Print details as each file is processed")
    args = parser.parse_args()
    csv2sql(args.directory, outfile=args.outfile, verbose=args.verbose)
