from astropy.coordinates import Angle
import astropy.units as u
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
