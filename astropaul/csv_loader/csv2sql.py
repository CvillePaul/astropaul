import argparse
import configparser
import inspect
from pathlib import Path
import sys

from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
import networkx as nx
import pandas as pd
import sqlalchemy as sa


class DataTransformation:
    def __init__(self, example_data: pd.DataFrame):
        self.example_data = example_data

    def get_sql_columns(self) -> list[tuple[str, object]]:
        return [(column_name, sa.String) for column_name in self.example_data.columns]

    def do_transformation(self, source: pd.DataFrame, dest: pd.DataFrame) -> pd.DataFrame:
        dest = pd.concat([dest, source[self.example_data.columns]], axis=1)
        return dest


class InferredTypeTransformation(DataTransformation):

    def __init__(self, example_data: pd.DataFrame):
        super().__init__(example_data)
        self.column_types = {}
        for column_name in self.example_data.columns:
            column_type = InferredTypeTransformation._series_to_column_type(self.example_data[column_name])
            self.column_types[column_name] = column_type

    @staticmethod
    def _series_to_column_type(series: pd.Series) -> str:
        if pd.api.types.is_integer_dtype(series):
            return sa.Integer
        elif pd.api.types.is_float_dtype(series):
            return sa.Float
        else:
            return sa.String  # Default to string for categorical/text data

    def get_sql_columns(self):
        return [(column_name, column_type) for column_name, column_type in self.column_types.items()]

    def do_transformation(self, source, dest):
        # verify that the data type of each column matches the example data provided during initialization
        for column_name, expected_column_type in self.column_types.items():
            source_column_type = InferredTypeTransformation._series_to_column_type(source[column_name])
            if source_column_type != expected_column_type:
                raise ValueError(f"Column {column_name} should be type {expected_column_type} but seeing {source_column_type}")
        return super().do_transformation(source, dest)


class CoordinateDataTransformation(DataTransformation):

    def __init__(self, example_data: pd.DataFrame):
        super().__init__(example_data)
        if len(self.example_data.columns) != 2:
            raise ValueError("Must pass in exactly two columns (for RA & Dec) for this transformation")
        self.ra_decimal, self.dec_decimal = self.example_data.columns
        self.ra_sexagesimal = self.ra_decimal + "_hms"
        self.dec_sexagesimal = self.dec_decimal + "_dms"

    @staticmethod
    def _get_angle_from_string(input: list[str], unit) -> Angle:
        angles = []
        for item in input:
            if ":" in item:
                angle = Angle(item, unit=unit)
            else:
                angle = Angle(item, unit=u.deg)
            angles.append(angle)
        return Angle(angles)

    def get_sql_columns(self):
        return [
            (self.ra_decimal, sa.Float),
            (self.dec_decimal, sa.Float),
            (self.ra_sexagesimal, sa.String),
            (self.dec_sexagesimal, sa.String),
        ]

    def do_transformation(self, source, dest):
        ra_angle = CoordinateDataTransformation._get_angle_from_string(list(source[self.ra_decimal]), unit=u.hourangle)
        dec_angle = CoordinateDataTransformation._get_angle_from_string(source[self.dec_decimal], unit=u.deg)
        dest[self.ra_decimal] = ra_angle.deg
        dest[self.dec_decimal] = dec_angle.deg
        dest[self.ra_sexagesimal] = ra_angle.to_string(decimal=False, sep=":", pad=True)
        dest[self.dec_sexagesimal] = dec_angle.to_string(decimal=False, sep=":", alwayssign=True, pad=True)
        return dest


class TableConfig:
    def __init__(self, config_file: Path, example_data: pd.DataFrame):
        self.table_name = string_to_db_style(config_file.name[:-4])
        self.has_id_column = False
        self.transformations = []
        self.constraints = {}  # key is dependent column, value is foreign (table, column)
        self.constraint_options = set()
        self.units = {}  # key is column name, value is string name of a class in astropy.units
        specified_transformations = {}  # key is frozen set of input columns, value is transformation class
        if config_file.exists():
            config = configparser.ConfigParser()
            config.optionxform = str  # don't convert keys to lower case
            config.read(config_file)
            # process table-wide options
            if specified_name := config.get("options", "table name", fallback=None):
                self.table_name = string_to_db_style(specified_name)
            self.has_id_column = config.getboolean("options", "create id column", fallback=False)
            all_constraint_options = {"log", "skip"}
            constraint_policy = config.get("options", "constraint policy", fallback=None)
            if constraint_policy:
                self.constraint_options = set([option.lower().strip() for option in constraint_policy.split(",")])
                if len(self.constraint_options) > 0 and not self.constraint_options.issubset(all_constraint_options):
                    raise ValueError(f"Constraint policy in {config_file} must be subset of: {all_constraint_options}")
            # process things that determine column order, transformations, type
            unspecified_columns = set(example_data.columns)

            if config.has_section("units"):
                self.units = {string_to_db_style(column_name): unit for column_name, unit in config.items("units")}
            if config.has_section("columns"):
                for column_name, column_type in config.items("columns"):
                    column_name = string_to_db_style(column_name)
                    if not column_name in unspecified_columns:
                        raise ValueError(f"Column type specified for unknown column {column_name} in {config_file}")
                    if column_type != "str":
                        raise ValueError(f"Unknown column type {column_type} for column {column_name} in file {config_file}")
                    specified_transformations[tuple([column_name])] = DataTransformation
                    unspecified_columns.remove(column_name)
            if config.has_section("transformations"):
                current_module = sys.modules[__name__]
                for transformation_name, column_list in config.items("transformations"):
                    transformation_class = getattr(current_module, transformation_name, None)
                    if not (
                        transformation_class
                        and inspect.isclass(transformation_class)
                        and transformation_class.__module__ == __name__
                    ):
                        raise ValueError(f"Unknown column transformer {transformation_name}")
                    transformation_columns = tuple([column_name.lower().strip() for column_name in column_list.split(",")])
                    specified_transformations[transformation_columns] = transformation_class
                    unspecified_columns.difference_update(transformation_columns)
            # now handle config items dealing with data relations
            if config.has_section("constraints"):
                for dependent_column, foreign_target in config.items("constraints"):
                    parts = foreign_target.split(".")
                    if not parts or len(parts) != 2:
                        raise ValueError(f"Bad constraint {foreign_target} for column {dependent_column} in {config_file}")
                    column_name = string_to_db_style(dependent_column)
                    if not column_name in example_data.columns:
                        raise ValueError(f"Constraint specified for unknown column {dependent_column} in {config_file}")
                    foreign_table = string_to_db_style(parts[0])
                    foreign_column = string_to_db_style(parts[1])
                    self.constraints[column_name] = (foreign_table, foreign_column)
        processed_columns = set()
        for column_name in example_data.columns:
            if column_name in processed_columns:
                continue  # this column is already handled by an existing transformation
            column_was_specified = False
            for affected_columns, transformation_class in specified_transformations.items():
                if column_name in affected_columns:
                    transformation = transformation_class(example_data[list(affected_columns)])
                    column_was_specified = True
                    processed_columns.update(affected_columns)
            if not column_was_specified:
                transformation = InferredTypeTransformation(example_data[[column_name]])
                processed_columns.add(column_name)
            self.transformations.append(transformation)
        pass


def string_to_db_style(name: str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return name.lower().replace(" ", "_")


def db_style_to_string(name: str) -> str:
    """Converts a name in database style back to a more pretty, human form"""
    answer = name.replace("_", " ")
    answer = answer.title()
    for string in ["Jd", "Utc", "Id", "Dssi", "Ra", "Hms", "Dms", "Pepsi", "Rv", "Pm"]:
        answer = answer.replace(string, string.upper())
    return answer


def create_table(metadata: sa.MetaData, table_config: TableConfig, table_metadata: pd.DataFrame) -> int:
    table_columns = []
    if table_config.has_id_column:
        table_columns.append(sa.Column("id", sa.Integer, primary_key=True, autoincrement=True))
    for transformation in table_config.transformations:
        for sql_name, sql_type in transformation.get_sql_columns():
            table_columns.append(sa.Column(sql_name, sql_type))
    sa.Table(table_config.table_name, metadata, *table_columns)
    for column_name, unit in table_config.units.items():
        table_metadata.loc[len(table_metadata)] = [table_config.table_name, column_name, "unit", unit]
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
            message = f"In file {file_name}, column {child_column} contains invalid values {constraint_violations}"
            if "log" in table_config.constraint_options:
                print(message)
            if "skip" in table_config.constraint_options:
                violation_mask = df_from_csv[child_column].isin(constraint_violations)
                df_from_csv = df_from_csv[~violation_mask]
            else:
                raise ValueError(message)
        # now build the DataFrame we'll convert to sql
        df_for_sql = pd.DataFrame()
        for transformation in table_config.transformations:
            df_for_sql = transformation.do_transformation(df_from_csv, df_for_sql)
        df_for_sql.to_sql(table_config.table_name, engine, if_exists="append", index=False)
        num_rows += len(df_from_csv)
    return num_rows


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


def csv2sql(base_dir: str, outfile: str, verbose: bool = False):
    engine = None
    try:
        outfile_path = Path(outfile)
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
            config_file = table_dir / f"{table_dir.name}.ini"
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
                print(f"Table {table_config.table_name}: ", end="")
            data_files = data_files_by_table[table_name]
            num_rows = insert_csv_data(engine, metadata, data_files, table_config)
            if verbose:
                print(f"{num_rows} rows inserted from {len(data_files)} files")

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
