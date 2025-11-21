import pandas as pd
import sqlalchemy as sa

from .db_utils import string_to_db_style

class KeyGenerator:

    def get_key_type(self):
        pass

    def generate_primary_keys(self, data: pd.DataFrame):
        pass


class SequentialIntKeyGenerator(KeyGenerator):
    def __init__(self):
        self.id_num = 0

    def get_key_type(self):
        return sa.Integer

    def generate_primary_keys(self, data: pd.DataFrame):
        num_rows = len(data)
        keys = [int(num) for num in range(self.id_num, self.id_num + num_rows)]
        self.id_num += num_rows
        return keys


class ColumnBasedKeyGenerator(KeyGenerator):
    def __init__(self, config_string: str):
        # TODO: parse config string & save details
        self.column_names = [string_to_db_style(column_name) for column_name in config_string.split(",")]

    def get_key_type(self):
        return sa.String

    def generate_primary_keys(self, data: pd.DataFrame):
        # TODO: implement keys from multiple columns and column substrings instead of single column
        return data[self.column_names].apply(lambda x: " ".join(x), axis=1)


class PepsiObservationKeyGenerator(KeyGenerator):
    def __init__(self, config_string: str):
        # TODO: parse config string & save details
        self.column_names = [string_to_db_style(column_name) for column_name in config_string.split(",")]

    def get_key_type(self):
        return sa.String

    def generate_primary_keys(self, data: pd.DataFrame):
        # TODO: implement keys from multiple columns and column substrings instead of single column
        return data[self.column_names].apply(lambda x: " ".join(x), axis=1)
