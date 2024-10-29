import pandas as pd
import astroplan as ap
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
import astropy.units as u
import numpy as np

class TargetListCreator:
    def __init__(self, **kwargs):
        self.kwargs = kwargs.copy()
        if steps := kwargs.get("steps"):
            self.steps = steps
            del self.kwargs["steps"]
        else:
            self.steps = []

    def calculate(self, initial_df: pd.DataFrame = None, steps=None, **kwargs) -> pd.DataFrame:
        """
        Create a new target list by running through self.steps and returning the result
        """
        intermediate_df = pd.DataFrame() if initial_df is None else initial_df
        merged_kwargs = {**self.kwargs, **kwargs}
        if steps is None:
            steps = self.steps
        for step in steps:
            intermediate_df = step(intermediate_df, **merged_kwargs)
        return intermediate_df


# following are methods to be used as steps of a TargetListCreator


def add_targets(df: pd.DataFrame, column_prefix="target_", **kwargs) -> pd.DataFrame:
    """
    Add rows for each target in the database
    this function intended always to be the first on in a chain: incoming df parameter ignored
    """
    targets = pd.read_sql(
        "select * from tom_target t;",
        kwargs["connection"],
        index_col="id",
    )
    new_column_names = {column: f"{column_prefix}{column}" for column in targets.columns}
    targets.rename(columns=new_column_names, inplace=True)
    return targets


def add_speckle(df: pd.DataFrame, column_prefix="speckle_", **kwargs) -> pd.DataFrame:
    count_column = f"{column_prefix}count"
    speckle = pd.read_sql("select * from tom_specklerawdata;", kwargs["connection"], index_col="target_id")
    speckle = speckle.groupby("target_id").id.agg(foo="count")
    speckle.rename(columns={"foo": count_column}, inplace=True)
    df = df.join(speckle)
    df[count_column] = df[count_column].fillna(0)
    df[f"{column_prefix}count"] = df[f"{column_prefix}count"].astype(int)
    return df


def add_spectra(df: pd.DataFrame, column_prefix="pepsi_", **kwargs) -> pd.DataFrame:
    count_column = f"{column_prefix}count"
    pepsi = pd.read_sql(
        "select * from tom_spectrumrawdata;",
        kwargs["connection"],
        index_col="target_id",
    )
    pepsi = pepsi.groupby("target_id").id.agg(foo="count")
    pepsi.rename(columns={"foo": count_column}, inplace=True)
    df = df.join(pepsi)
    df[count_column] = df[count_column].fillna(0)
    df[f"{column_prefix}count"] = df[f"{column_prefix}count"].astype(int)
    return df


def add_lists(df: pd.DataFrame, column_prefix="list_", **kwargs) -> pd.DataFrame:
    conn = kwargs["connection"]
    for (target_list,) in conn.execute("select name from tom_targetlist;").fetchall():
        list_members = [
            result[0]
            for result in conn.execute(
                """
                select t.name
                from tom_targetlist tl
                join tom_targetlist_targets tlt on tlt.targetlist_id = tl.id
                join tom_target t on t.id = tlt.target_id
                where tl.name = ?
                ;
                """,
                [target_list],
            ).fetchall()
        ]
        column_name = f'{column_prefix}{target_list.replace(" ", "_")}'
        df[column_name] = df.target_name.isin(list_members)
    return df


def add_ephemerides(df: pd.DataFrame, column_prefix="ephem_", **kwargs) -> pd.DataFrame:
    ephem = pd.read_sql(
        """
        select t.id, bp.system, bp.member, bp.period, bp.t0, bp.duration, bp.depth
        from tom_binaryparameters bp
        join tom_scienceresult sr on sr.id = bp.scienceresult_ptr_id
        join tom_target t on t.id = sr.target_id
        ;""",
        kwargs["connection"],
        index_col="id",
    )
    ephem["period"] = ephem["period"] / 3600 / 24 # convert from seconds to days for convenience
    ephem["duration"] = ephem["duration"] / 3600 #/ 24 # convert from seconds to hours for convenience
    new_column_names = {column: f"{column_prefix}{column}" for column in ephem.columns}
    ephem.rename(columns=new_column_names, inplace=True)
    df = df.join(ephem)
    return df


def add_tess(df: pd.DataFrame, column_prefix="tess_", **kwargs) -> pd.DataFrame:
    tess = pd.read_sql(
        """
        select ca.target_id, tt.*
        from tom_tess_ticv8 tt
        join tom_catalogassociation ca on ca.catalog_id = tt.Identifier
        where ca.association = 'Primary ID' and ca.catalog = 'TESS TICv8'
        """,
        kwargs["connection"],
        index_col="target_id",
    )
    tess.drop("id", axis=1, inplace=True)
    new_column_names = {column: f"{column_prefix}{column}" for column in tess.columns}
    tess.rename(columns=new_column_names, inplace=True)
    df = df.join(tess)
    return df


def add_coords(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Add fundamental things like coordinates from available catalog data"""
    tess_mapping = {
        "ra": "tess_ra",
        "dec": "tess_dec",
        "pmra": "tess_pmRA",
        "pmdec": "tess_pmDEC",
        "parallax": "tess_plx",
        "Vmag": "tess_Vmag",
        "Teff": "tess_Teff",
    }
    if "tess_version" in df.columns:
        for field, tess_field in tess_mapping.items():
            df[field] = df[tess_field]
    # add sexagesimal versions of the coordinates for convenience
    df["RA"] = [Angle(ra, unit=u.deg).to_string(unit=u.hour, decimal=False, precision=2, sep=":") for ra in df["ra"]]
    df["Dec"] = [
        Angle(dec, unit=u.deg).to_string(unit=u.deg, decimal=False, precision=2, sep=":", alwayssign=True)
        for dec in df["dec"]
    ]
    return df


def hide_cols(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    cols_to_remove = []
    if prefix := kwargs.get("prefix"):
        prefix_cols = [col for col in df.columns if col.startswith(prefix)]
        cols_to_remove += prefix_cols
    if len(cols_to_remove) > 0:
        df.drop(labels=cols_to_remove, axis=1, inplace=True)
    return df


def add_observability(
    df: pd.DataFrame,
    observer: ap.Observer,
    time_segment: tuple[Time, Time], # if more than a day apart, calculates observing nights during given range 
    constraints: list[ap.Constraint] = None,
    column_prefix: str = "",
    **kwargs,
) -> pd.DataFrame:
    if constraints is None:
        constraints = [
            ap.AltitudeConstraint(30 * u.deg, 80 * u.deg),
            ap.AirmassConstraint(2),
            # ap.AtNightConstraint.twilight_civil(),
        ]
    targets = ap.FixedTarget(coord=SkyCoord(ra=df["ra"].values * u.deg, dec=df["dec"].values * u.deg))
    # astroplan needs a contiguous range of time to test observability
    # break given range into observing nights & test separately, combining results
    beg_time, end_time = time_segment
    sample_intervals = [] # list of tuples of time range(s) to test
    if (end_time - beg_time) > 1 * u.day:
        curr_time = beg_time
        term_time = observer.sun_rise_time(end_time, which="next")
        while curr_time < term_time:
            beg_night = observer.sun_set_time(curr_time, which="nearest")
            end_night = observer.sun_rise_time(beg_night, which="next")
            sample_intervals.append((beg_night, end_night))
            curr_time += 1 * u.day
    else:
        sample_intervals = [time_segment]
    
    observability = np.array([False] * len(df))
    for sample_interval in sample_intervals:
        print(sample_interval[0].iso, sample_interval[1].iso)
        interval_observability = np.array(ap.is_observable(constraints, observer, targets, sample_interval))
        observability |= interval_observability
    df[f"{column_prefix}observable"] = observability
    return df


# def add_observability(df: pd.DataFrame, column_prefix="observable_", session_beg: Time=None, session_end: Time=None, **kwargs) -> pd.DataFrame:
#     # make astroplan objects for each target
#     fixed_targets = [
#         ap.FixedTarget(
#             coord=SkyCoord(
#                 frame="icrs",
#                 obstime=Time("2000.0", format="jyear", scale="tdb"),
#                 ra=row["ra"] * u.deg,
#                 dec=row["dec"] * u.deg,
#                 pm_ra_cosdec=row["pmra"] * u.mas / u.yr,
#                 pm_dec=row["pmdec"] * u.mas / u.yr,
#                 # distance=row["parallax"] * u.pc,
#             ),
#             name=row["target_name"],
#         )
#         for (_, row) in df.iterrows()
#     ]

#     constraints = [
#         ap.AltitudeConstraint(30 * u.deg, 80 * u.deg),
#         ap.AirmassConstraint(2),
#         # ap.AtNightConstraint.twilight_civil(),
#     ]

#     # calculate observability
#     df[f"{column_prefix}observable"] = ap.is_observable(
#         constraints, observer, fixed_targets, (session_beg, session_end)
#     )

#     # calculate rise/meridian/set times
#     with warnings.catch_warnings():
#         # ignore warnings thrown when targets are never or always visible
#         warnings.simplefilter("ignore", AstropyWarning)
#         rise_times = observer.target_rise_time(
#             session_beg, fixed_targets, which="nearest"
#         ).filled(session_beg)
#         df[f"{column_prefix}_rise_jd"] = rise_times.jd
#         df[f"{column_prefix}_rise_time"] = list(
#             map(
#                 lambda x: x.to_datetime(timezone=observer.timezone).isoformat(
#                     sep=" ", timespec="seconds"
#                 ),
#                 rise_times,
#             )
#         )
#         meridian_times = observer.target_meridian_transit_time(
#             session_beg, fixed_targets, which="nearest"
#         )
#         df[f"{column_prefix}meridian_jd"] = meridian_times.jd
#         df[f"{column_prefix}meridian_time"] = list(
#             map(
#                 lambda x: x.to_datetime(timezone=observer.timezone).isoformat(
#                     sep=" ", timespec="seconds"
#                 ),
#                 meridian_times,
#             )
#         )
#         set_times = observer.target_set_time(
#             session_beg, fixed_targets, which="nearest"
#         ).filled(session_end)
#         df[f"{column_prefix}set_jd"] = set_times.jd
#         df[f"{column_prefix}set_time"] = list(
#             map(
#                 lambda x: x.to_datetime(timezone=observer.timezone).isoformat(
#                     sep=" ", timespec="seconds"
#                 ),
#                 set_times,
#             )
#         )
#         return df

# # TODO: add columns for beginning, ending, and max airmass or altitude
