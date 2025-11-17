import configparser
import inspect
import sys

import pandas as pd

from .DataTransformation import *
from .db_utils import *


class TableConfig:
    def __init__(self, config_file: Path, example_data: pd.DataFrame):
        self.table_name = string_to_db_style(config_file.name[:-4])
        self.has_id_column = False
        self.transformations = []
        self.constraints = {}  # key is dependent column, value is foreign (table, column)
        self.constraint_options = set()
        self.units = {}  # key is column name, value is a valid astropy.units string
        self.metadata = None
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
            if self.table_name:
                pass
            if config.has_section("metadata"):
                self.metadata = {key: value for key, value in config.items("metadata")}
            if config.has_section("units"):
                self.units = {string_to_db_style(column_name): unit for column_name, unit in config.items("units")}

            # process things that determine column order, transformations, type
            unspecified_columns = set(example_data.columns)
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
