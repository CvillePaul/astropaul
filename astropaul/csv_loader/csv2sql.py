import argparse
import configparser
from functools import partial
from pathlib import Path
import typing

from astropy.coordinates import Angle
import astropy.units as u
import networkx as nx
import pandas as pd
import sqlalchemy as sa

ColumnTransformer = typing.Callable[[str, pd.DataFrame, pd.DataFrame], pd.DataFrame]


def DefaultColumnTransformer(column_name: str, source: pd.DataFrame, dest: pd.DataFrame) -> pd.DataFrame:
    dest_column_name = db_style_to_string(column_name)
    dest[dest_column_name] = source[column_name]
    return dest


def SexagesimalOrDecimalDegreesColumnTransformer(
    column_name: str, source: pd.DataFrame, dest: pd.DataFrame, hourangle: bool
) -> pd.DataFrame:
    if source.empty:
        raise ValueError("Source DataFrame is empty")
    dest_column_name = db_style_to_string(column_name)
    sample_value = source.iloc[0][column_name]
    if ":" in sample_value:
        unit = u.hourangle if hourangle else u.deg
        data = Angle(source[column_name], unit=unit)
    else:
        data = Angle(source[column_name], unit=u.deg)
    dest[dest_column_name] = data.deg
    dest[f"{dest_column_name} (HMS)"] = data.to_string(decimal=False, sep=":")


class TableConfig:
    def __init__(self):
        self.table_name: str = "table"
        self.has_id_column: bool = False
        self.constraints: dict[str, tuple[str, str]] = dict()  # column name: (table name, column name) of foreign table
        self.column_types: dict[str, object] = dict()  # manual column type overrides
        self.column_transformations: dict[str, ColumnTransformer] = (
            dict()
        )  # how to change the csv contents into some other form

    def validate_csv_schemas(self, data_files: list[tuple[Path, pd.DataFrame]]):
        reference_file, reference_df = data_files[0]  # Use first file as reference
        all_columns = set(reference_df.columns)
        specified_columns = {}
        for column, column_type in self.column_types.items():
            if column not in all_columns:
                raise ValueError(f"Column {column} specified in config file but not present in {reference_file}")
            specified_columns[column] = column_type
        inferred_columns = {
            column: series_to_column_type(reference_df[column])
            for column in reference_df.columns
            if column not in specified_columns
        }
        self.column_types = {**inferred_columns, **specified_columns}
        for column_name, column_type in self.column_types.items():
            if column_name in self.column_transformations:
                continue  # don't override manually specified column transformations
            self.column_transformations[column_name] = DefaultColumnTransformer
        for file_name, data_file in data_files[1:]:  # validate all remaining files against schema we just built
            df_columns = set(data_file.columns)
            if df_columns != all_columns:
                raise ValueError(
                    (
                        f"Schema mismatch in {file_name} vs reference file {reference_file}. "
                        "Expected columns: {all_columns}, Found: {df_columns}"
                    )
                )
            for column, column_type in inferred_columns.items():
                this_type = series_to_column_type(data_file[column])
                if this_type != column_type:
                    raise ValueError(
                        (f"Data type mismatch in {file_name} for column {column}: " "{this_type} instead of {column_type}")
                    )

    @classmethod
    def from_config_file(cls, config_file: Path) -> "TableConfig":
        answer = TableConfig()
        answer.table_name = config_file.name[:-4]
        if config_file.exists():
            config = configparser.ConfigParser()
            config.optionxform = str  # don't convert keys to lower case
            config.read(config_file)
            answer.table_name = config.get("options", "table name", fallback=answer.table_name)
            answer.has_id_column = config.getboolean("options", "create id column", fallback=False)
            if config.has_section("columns"):
                for column_name, column_type_str in config.items("columns"):
                    if not (column_type := string_to_column_type(column_type_str)):
                        raise ValueError(f"Unknown column type {column_type_str} in file {config_file}")
                    answer.column_types[column_name] = column_type
            if config.has_section("constraints"):
                for column_name, foreign_column in config.items("constraints"):
                    parts = foreign_column.split(".")
                    if not parts or len(parts) != 2:
                        raise ValueError(f"Bad constraint {foreign_column} for column {column_name} in {config_file}")
                    foreign_table, foreign_column = parts
                    answer.constraints[column_name] = (foreign_table, foreign_column)
            if config.has_section("transformations"):
                for column_name, transformer_name in config.items("transformations"):
                    match transformer_name.lower():
                        case "hms_or_decimal_degrees":
                            transformer = partial(SexagesimalOrDecimalDegreesColumnTransformer, hourangle=True)
                        case "dms_or_decimal_degrees":
                            transformer = partial(SexagesimalOrDecimalDegreesColumnTransformer, hourangle=False)
                        case _:
                            raise ValueError(
                                f"Unknown column transformation {transformer_name} for column {column_name} of {answer.table_name}"
                            )
                    answer.column_transformations[column_name] = transformer
        return answer


def string_to_db_style(name: str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return name.lower().replace(" ", "_")

def db_style_to_string(name: str) -> str:
    """Converts a name in database style back to a more pretty, human form"""
    answer = name.replace("_", " ")
    answer = answer.title()
    for string in ["Jd", "Utc", "Id", "Dssi", "Ra", "Hms", "Dms"]:
        answer = answer.replace(string, string.upper())
    return answer

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


def create_table(metadata, table_config: TableConfig):
    table_columns = []
    if table_config.has_id_column:
        table_columns.append(sa.Column("id", sa.Integer, primary_key=True, autoincrement=True))
    for column_name, column_type in table_config.column_types.items():
        table_columns.append(sa.Column(string_to_db_style(column_name), column_type))
    table = sa.Table(string_to_db_style(table_config.table_name), metadata, *table_columns)
    return table


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
        source_table = metadata.tables[string_to_db_style(source_table_name)]
        with engine.connect() as conn:
            result = conn.execute(sa.select(source_table.c[string_to_db_style(source_column)]))
            validation_values[child_column] = set([row[0] for row in result])
    # validate and insert rows from each data file
    num_rows = 0
    for file_name, data_file in data_files:
        for child_column in table_config.constraints.keys():
            valid_column_values = validation_values[child_column]
            for value in data_file[child_column].unique():
                if not value in valid_column_values:
                    raise ValueError(
                        f"In file {file_name} for table {table_config.table_name}, value {value} not found in {source_table_name}.{source_column}"
                    )
        db_frame = data_file.copy()
        db_frame.columns = [string_to_db_style(column) for column in db_frame.columns]
        db_frame.to_sql(string_to_db_style(table_config.table_name), engine, if_exists="append", index=False)
        num_rows += len(data_file)
    return num_rows


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
        table_configs = {}  # key is table name, value is TableConfig
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
            config_file = table_dir / f"{table_dir.name}.ini"
            table_config = TableConfig.from_config_file(config_file)
            table_config.validate_csv_schemas(data_files)
            table_configs[table_config.table_name] = table_config
            create_table(metadata, table_config)
            if verbose:
                print(f"Created table {table_name} with {len(table_config.column_types)} columns")
        metadata.create_all(engine)

        # load data into the new tables
        table_order = determine_table_order(table_configs)
        for table_name in table_order:
            data_files = data_files_by_table[table_name]
            table_config = table_configs[table_name]
            num_rows = insert_csv_data(engine, metadata, data_files, table_config)
            if verbose:
                print(f"{num_rows:5d} rows inserted into table {table_config.table_name} from {len(data_files)} files")

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
