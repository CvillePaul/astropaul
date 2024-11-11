from dataclasses import dataclass

import pandas as pd
import astroplan as ap
from astropy.coordinates import SkyCoord, Angle
from astropy.time import Time
import astropy.units as u
import numpy as np


class TargetList:
    def __init__(
        self,
        target_list: pd.DataFrame = None,
        column_groups: dict[str, tuple[list[str], list[str]]] = None,
        other_lists: dict[str, pd.DataFrame] = None,
    ):
        self.target_list = target_list if target_list is not None else pd.DataFrame()
        self.colum_groups = column_groups or {}
        self.other_lists = other_lists or {}

    def copy(self) -> "TargetList":
        answer = TargetList(
            self.target_list.copy(),
            self.colum_groups.copy(),
            {key: value.copy() for key, value in self.other_lists.items()},
        )
        return answer

    @classmethod
    def merge(
        cls,
        first: "TargetList",
        target_list: pd.DataFrame = pd.DataFrame(),
        column_groups: dict[str, tuple[list[str], list[str]]] = {},
        other_lists: dict[str, pd.DataFrame] = {},
    ) -> "TargetList":
        first = first if first else TargetList()
        answer = TargetList(first.target_list, first.colum_groups, first.other_lists)
        if answer.target_list.empty or answer.target_list is None:
            answer.target_list = target_list
        else:
            answer.target_list = answer.target_list.join(target_list)
        answer.colum_groups = {**answer.colum_groups, **column_groups}
        for name, other_list in other_lists.items():
            answer.other_lists[name] = other_list
        return answer


class TargetListCreator:
    def __init__(self, **kwargs):
        self.kwargs = kwargs.copy()
        if steps := kwargs.get("steps"):
            self.steps = steps
            del self.kwargs["steps"]
        else:
            self.steps = []

    def calculate(self, initial_df: TargetList = TargetList(), steps=None, **kwargs) -> TargetList:
        """
        Create a new target list by running through self.steps and returning the result
        """
        intermediate_tl = initial_df if initial_df else TargetList()
        merged_kwargs = {**self.kwargs, **kwargs}
        if steps is None:
            steps = self.steps
        for step in steps:
            intermediate_tl = step(intermediate_tl, **merged_kwargs)
        return intermediate_tl


# following are methods to be used as steps of a TargetListCreator


def add_targets(tl: TargetList = TargetList(), column_prefix="Target ", **kwargs) -> TargetList:
    """
    Add rows for each target in the database
    """
    targets = pd.read_sql(
        "select * from tom_target t;",
        kwargs["connection"],
        index_col="id",
    )
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in targets.columns}
    targets.rename(columns=new_column_names, inplace=True)
    column_groups = {"Target": (["Name"], ["Source", "Type"])}
    return TargetList.merge(tl, targets, column_groups, {})


def add_speckle(tl: TargetList, column_prefix="Speckle ", **kwargs) -> pd.DataFrame:
    # first, add a column for total observations to main table
    count_column = f"{column_prefix}Count"
    speckle_count = pd.read_sql(
        f"""
        select target_id, count(id) as '{count_column}'
        from tom_specklerawdata
        group by target_id
        ;
        """,
        kwargs["connection"],
        index_col="target_id",
    )
    column_groups = {"Speckle Count": ([count_column], [])}
    # next, add a separate table of all speckle observations
    speckle = pd.read_sql("select * from tom_specklerawdata", kwargs["connection"], index_col="target_id")
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in speckle.columns}
    speckle.rename(columns=new_column_names, inplace=True)
    answer = TargetList.merge(tl, speckle_count, column_groups, {"Speckle": speckle})
    answer.target_list[count_column] = answer.target_list[count_column].fillna(int(0))
    return answer


def add_spectra(tl: TargetList, column_prefix="PEPSI ", **kwargs) -> TargetList:
    # first, add a column for total observations to main table
    count_column = f"{column_prefix}Count"
    spectra_count = pd.read_sql(
        f"""
        select target_id, count(id) as '{count_column}'
        from tom_spectrumrawdata
        group by target_id
        ;""",
        kwargs["connection"],
        index_col="target_id",
    )
    column_groups = {"Spectra Count": ([count_column], [])}
    # next, add a separate table of all spectra
    spectra = pd.read_sql("select * from tom_spectrumrawdata", kwargs["connection"], index_col="target_id")
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in spectra.columns}
    spectra.rename(columns=new_column_names, inplace=True)
    answer = TargetList.merge(tl, spectra_count, column_groups, {"Spectra": spectra})
    answer.target_list[count_column] = answer.target_list[count_column].fillna(int(0))
    return answer


def add_lists(tl: TargetList, column_prefix="List ", **kwargs) -> TargetList:
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
                ;""",
                [target_list],
            ).fetchall()
        ]
        column_name = f'{column_prefix}{target_list.replace(" ", "_").capitalize()}'
        answer = tl.copy()
        # TODO: remove hard-coded column name in following line
        answer.target_list[column_name] = answer.target_list["Target Name"].isin(list_members)
    return answer


def add_ephemerides(tl: TargetList, column_prefix="Ephem ", **kwargs) -> TargetList:
    ephem = pd.read_sql(
        """
        select t.id, t.name, bp.system, bp.member, bp.period, bp.t0, bp.duration, bp.depth
        from tom_binaryparameters bp
        join tom_scienceresult sr on sr.id = bp.scienceresult_ptr_id
        join tom_target t on t.id = sr.target_id
        ;""",
        kwargs["connection"],
        index_col="id",
    )
    ephem["period"] = ephem["period"] / 3600 / 24  # convert from seconds to days for convenience
    ephem["duration"] = ephem["duration"] / 3600  # convert from seconds to hours for convenience
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in ephem.columns}
    ephem.rename(columns=new_column_names, inplace=True)
    answer = tl.copy()
    answer.other_lists[column_prefix.strip().replace("_", "")] = ephem
    return answer


def add_tess(tl: TargetList, column_prefix="TESS ", **kwargs) -> TargetList:
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
    identifier_column = f"{column_prefix}Identifier"
    columns = list(new_column_names.values())
    columns.remove(identifier_column)
    column_groups = {column_prefix: ([identifier_column], columns)}
    return TargetList.merge(tl, tess, column_groups, {})


def add_coords(tl: TargetList, **kwargs) -> TargetList:
    """Add fundamental things like coordinates if available from catalog data"""
    # TODO: remove hard coded column prefix from this dictionary
    tess_mapping = {
        "ra": "TESS ra",
        "dec": "TESS dec",
        "pmra": "TESS pmRA",
        "pmdec": "TESS pmDEC",
        "parallax": "TESS plx",
        "Vmag": "TESS Vmag",
        "Teff": "TESS Teff",
    }
    answer = tl.copy()
    if "TESS version" in answer.target_list.columns:
        for field, tess_field in tess_mapping.items():
            answer.target_list[field] = answer.target_list[tess_field]
    # add sexagesimal versions of the coordinates for convenience
    answer.target_list["RA"] = [
        Angle(ra, unit=u.deg).to_string(unit=u.hour, decimal=False, precision=2, sep=":") for ra in answer.target_list["ra"]
    ]
    answer.target_list["Dec"] = [
        Angle(dec, unit=u.deg).to_string(unit=u.deg, decimal=False, precision=2, sep=":", alwayssign=True)
        for dec in answer.target_list["dec"]
    ]
    answer.colum_groups = {
        **answer.colum_groups,
        "Coordinates": (["RA", "Dec", "ra", "dec", "Vmag", "Teff"], ["pmra", "pmdec", "parallax"]),
    }
    return answer


def hide_cols(tl: TargetList, **kwargs) -> TargetList:
    answer = tl.copy()
    if prefix := kwargs.get("prefix"):
        cols_to_remove = [col for col in tl.target_list.columns if col.startswith(prefix)]
        if len(cols_to_remove) > 0:
            answer.target_list.drop(labels=cols_to_remove, axis=1, inplace=True)
    del answer.colum_groups[prefix]
    return answer


def add_observability(
    tl: TargetList,
    observer: ap.Observer,
    time_segment: tuple[Time, Time],  # if more than a day apart, calculates observing nights during given range
    constraints: list[ap.Constraint] = None,
    column_prefix: str = "Observable ",
    calc_max_altitude: bool = False,
    **kwargs,
) -> TargetList:
    if constraints is None:
        constraints = [
            ap.AltitudeConstraint(30 * u.deg, 80 * u.deg),
            ap.AirmassConstraint(2),
        ]
    coords = SkyCoord(ra=tl.target_list["ra"].values * u.deg, dec=tl.target_list["dec"].values * u.deg)
    targets = ap.FixedTarget(coord=coords)
    # astroplan needs a contiguous range of time to test observability
    # break given range into observing nights & test separately, combining results
    beg_time, end_time = time_segment
    sample_intervals = []  # list of tuples of time range(s) to test
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

    # calculate observability for each nightly time range
    answer = tl.copy()
    overall_observability = np.array([False] * len(answer.target_list))
    overall_max_alts = np.array([-90.0] * len(answer.target_list))
    nightly_columns = []
    for sample_interval in sample_intervals:
        interval_observability = np.array(ap.is_observable(constraints, observer, targets, sample_interval))
        column_name = f"{column_prefix}{sample_interval[0].iso[:10]}"
        answer.target_list[column_name] = interval_observability
        overall_observability |= interval_observability
        nightly_columns.append(column_name)
        if calc_max_altitude:
            column_name += " Max Alt"
            time_grid = ap.utils.time_grid_from_range(sample_interval, time_resolution=15 * u.min)
            alts = [np.max(observer.altaz(time_grid, coord).alt.value) for coord in coords]
            answer.target_list[column_name] = alts
            nightly_columns.append(column_name)
            overall_max_alts = [
                np.max([current_max_alt, this_alt]) for current_max_alt, this_alt in zip(overall_max_alts, alts)
            ]
    overall_column = f"{column_prefix}Any Night"
    answer.target_list[overall_column] = overall_observability
    main_columns = [overall_column]
    if calc_max_altitude:
        max_alt_column = f"{column_prefix}Max Alt"
        answer.target_list[max_alt_column] = overall_max_alts
        main_columns.append(max_alt_column)
    answer.colum_groups[column_prefix.strip()] = (main_columns, nightly_columns)
    return answer


def filter_targets(
    tl: TargetList,
    criteria=True,
    **kwargs,
) -> TargetList:
    answer = tl.copy()
    answer.target_list = answer.target_list[criteria]
    return answer
