import collections
import datetime
import platform
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, Angle, get_body
from astropy.time import Time
import astropy.units as u
import itables
import numpy as np
import pandas as pd


class ObservingSession:
    def __init__(
        self,
        observer: ap.Observer,
        observing_segments: list[tuple[Time, Time]] = [],
    ):
        self.observer = observer
        self.observing_segments = observing_segments

    def __repr__(self):
        if self.observer.name:
            answer = f"{self.observer.name}"
        else:
            answer = f"Site ({self.observer.latitude.value:.3f}, {self.observer.longitude.value:.3f})"
        answer += f" in {self.observer.timezone}\n"
        if self.observing_segments:
            for beg, end in self.observing_segments:
                answer += f"    {beg.iso[:19]} to {end.iso[:19]}\n"
        else:
            answer += "No observing segments defined\n"
        return answer

    def _determine_nighttime(self, night: Time) -> tuple[Time, Time]:
        beg_night = self.observer.sun_set_time(night, which="nearest")
        end_night = self.observer.sun_rise_time(beg_night, which="next")
        return (beg_night, end_night)

    def add_full_day(self, day: str | Time):
        self.observing_segments.append(self._determine_nighttime(Time(day)))

    def add_half_day(self, day: str | Time, first_half: bool = True):
        beg, end = self._determine_nighttime(day)
        mid = Time(beg.jd + ((end.jd - beg.jd) / 2), format="jd")
        if first_half:
            self.observing_segments.append((beg, mid))
        else:
            self.observing_segments.append((mid, end))

    def add_day_range(self, range: tuple[str, str] | tuple[Time, Time]):
        beg, end = Time(range[0]), Time(range[1])
        day = beg
        while day < end:
            self.add_full_day(day)
            day += 1 * u.day


class TargetList:
    def __init__(
        self,
        name: str = "Target List",
        target_list: pd.DataFrame = None,
        column_groups: dict[str, tuple[list[str], list[str]]] = None,
        other_lists: dict[str, pd.DataFrame] = {},
    ):
        self.name = name
        self.target_list = target_list if target_list is not None else pd.DataFrame()
        self.column_groups = column_groups or collections.defaultdict(list)
        self.other_lists = other_lists

    def copy(self) -> "TargetList":
        answer = TargetList(
            name=self.name,
            target_list=self.target_list.copy(),
            column_groups=self.column_groups.copy(),
            other_lists={key: value.copy() for key, value in self.other_lists.items()},
        )
        return answer

    def addColumns(self, columns: dict[str, Any] = {}):
        pass

    def addOther(self, name: str, other: pd.DataFrame):
        self.other_lists[name] = other

    def summarize(self) -> str:
        answer = f"{self.name}\n"
        target_types = collections.Counter(self.target_list["Target Type"])
        answer += f"{len(self.target_list)} targets:\n"
        for type, count in target_types.items():
            answer += f"    {count:4d} {type}\n"
        answer += "Columns (primary, secondary):\n"
        for name, (primary, secondary) in self.column_groups.items():
            answer += f"    {name}: ({len(primary)}, {len(secondary)})\n"
        answer += "Other tables:\n"
        for name, other in self.other_lists.items():
            answer += f"    {len(other):4d} {name}\n"
        return answer

    def to_html(self, file: str = None) -> str:
        buttons = [
            {
                "extend": "columnToggle",
                "columns": [f"{col}:title" for col in secondary_cols],
                "text": group,
                "redrawCalculations": True,
            }
            for group, (_, secondary_cols) in self.column_groups.items()
            if secondary_cols
        ]
        table_name = "targetList"
        caption = f"{self.name}, created {datetime.datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}"
        html = itables.to_html_datatable(
            self.target_list,
            caption=caption,
            table_id=table_name,
            connected=True,
            paging=False,
            showIndex=False,
            maxBytes=0,
            maxColumns=0,
            autoWidth=False,
            layout={"topStart": {"searchBuilder": True, "buttons": buttons}},
            classes="display nowrap compact cell-border",
            # columnDefs=[{"targets": ["Target Source:title"], "visible": False}],
        )
        html += f"""
            <style>
            #{table_name} th {{
                white-space: normal;
                word-wrap: break-word;
            }}
            caption {{
                text-align: left;
                font-style: italic;
            }}
            </style>
            """

        if file:
            with open(file, "w") as f:
                f.write(html)
        return html

    @classmethod
    def merge(
        cls,
        first: "TargetList",
        target_list: pd.DataFrame = pd.DataFrame(),
        column_groups: dict[str, tuple[list[str], list[str]]] = {},
        other_lists: dict[str, pd.DataFrame] = {},
    ) -> "TargetList":
        first = first if first else TargetList()
        answer = TargetList(
            name=first.name, target_list=first.target_list, column_groups=first.column_groups, other_lists=first.other_lists
        )
        if answer.target_list.empty or answer.target_list is None:
            answer.target_list = target_list
        else:
            answer.target_list = answer.target_list.join(target_list)
        answer.column_groups = {**answer.column_groups, **column_groups}
        for name, other_list in other_lists.items():
            answer.other_lists[name] = other_list
        return answer


class TargetListCreator:
    def __init__(self, name: str = "Standard", steps: list = None, **kwargs):
        self.name = name
        self.kwargs = kwargs.copy()
        self.steps = steps

    def calculate(self, initial_df: TargetList = None, steps=None, verbose: bool = False, **kwargs) -> TargetList:
        """
        Create a new target list by running through self.steps and returning the result
        """
        if not initial_df:
            initial_df = TargetList(name=self.name)
        intermediate_tl = initial_df if initial_df else TargetList(name=self.name)
        merged_kwargs = {**self.kwargs, **kwargs}
        if steps is None:
            steps = self.steps
        for step in steps:
            intermediate_tl = step(intermediate_tl, **merged_kwargs)
            if verbose:
                print(f"{len(intermediate_tl.target_list):4d} targets, {step=}")
        if verbose:
            print(intermediate_tl.summarize())
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
    column_groups = {"Target": ([f"{column_prefix}Name"], [f"{column_prefix}Source", f"{column_prefix}Type"])}
    return TargetList.merge(tl, targets, column_groups, {})


def concat_dataframe(tl: TargetList, other_df: pd.DataFrame, **kwargs) -> TargetList:
    return TargetList(
        name=tl.name,
        target_list=pd.concat([tl.target_list, other_df]),
        column_groups=tl.column_groups,
        other_lists=tl.other_lists,
    )


def add_speckle(tl: TargetList, column_prefix="Speckle ", **kwargs) -> TargetList:
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


def add_spectra(tl: TargetList, column_prefix="Spectra ", **kwargs) -> TargetList:
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
    list_names = [name[0] for name in conn.execute("select name from tom_targetlist;").fetchall()]
    answer = tl.copy()
    for target_list in list_names:
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
        column_name = f"{column_prefix}{target_list}"
        # column_name = f'{column_prefix}{target_list.replace(" ", "_").capitalize()}'
        # TODO: remove hard-coded column name in following line
        answer.target_list[column_name] = answer.target_list["Target Name"].isin(list_members)
        answer.column_groups[column_prefix] = ([], [f"{column_prefix}{name}" for name in list_names])
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
    answer.addOther(column_prefix.strip().replace("_", ""), ephem)
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
    objId_column = f"{column_prefix}objID"
    columns = list(new_column_names.values())
    columns.remove(objId_column)
    column_groups = {column_prefix: ([objId_column], columns)}
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
        Angle(ra, unit=u.deg).to_string(unit=u.hour, decimal=False, precision=2, sep=":", pad=True)
        for ra in answer.target_list["ra"]
    ]
    answer.target_list["Dec"] = [
        Angle(dec, unit=u.deg).to_string(unit=u.deg, decimal=False, precision=2, sep=":", alwayssign=True, pad=True)
        for dec in answer.target_list["dec"]
    ]
    answer.column_groups = {
        **answer.column_groups,
        "Coordinates": (["RA", "Dec", "ra", "dec", "Vmag", "Teff"], ["pmra", "pmdec", "parallax"]),
    }
    return answer


def hide_cols(tl: TargetList, **kwargs) -> TargetList:
    answer = tl.copy()
    if prefix := kwargs.get("prefix"):
        cols_to_remove = [col for col in tl.target_list.columns if col.startswith(prefix)]
        if len(cols_to_remove) > 0:
            answer.target_list.drop(labels=cols_to_remove, axis=1, inplace=True)
    del answer.column_groups[prefix]
    return answer


def add_observability(
    tl: TargetList,
    observing_session: ObservingSession,
    constraints: list[ap.Constraint] = None,
    column_prefix: str = "Observable ",
    calc_max_altitude: bool = False,
    calc_moon_distance: bool = False,
    **kwargs,
) -> TargetList:
    if not constraints:
        constraints = [
            ap.AltitudeConstraint(30 * u.deg, 80 * u.deg, True),
            # ap.AirmassConstraint(1),
        ]
    coords = SkyCoord(ra=tl.target_list["ra"].values * u.deg, dec=tl.target_list["dec"].values * u.deg)

    # calculate observability for each time segment of the observing session
    answer = tl.copy()
    overall_any_night = np.array([False] * len(answer.target_list))
    overall_every_night = np.array([True] * len(answer.target_list))
    overall_max_alts = np.array([-90.0] * len(answer.target_list))
    overall_min_moon_dist = np.array([180.0] * len(answer.target_list))
    nightly_columns = []
    for beg, end in observing_session.observing_segments:
        interval_observability = np.array(ap.is_observable(constraints, observing_session.observer, coords, time_range=(beg, end)))
        column_name = f"{column_prefix}{beg.iso[:10]}"
        answer.target_list[column_name] = interval_observability
        overall_any_night |= interval_observability
        overall_every_night &= interval_observability
        nightly_columns.append(column_name)
        if calc_max_altitude:
            max_alt_column = f"{column_name} Max Alt"
            time_grid = ap.utils.time_grid_from_range((beg, end), time_resolution=15 * u.min)
            dist = [np.max(observing_session.observer.altaz(time_grid, coord).alt.value) for coord in coords]
            answer.target_list[max_alt_column] = dist
            nightly_columns.append(max_alt_column)
            overall_max_alts = [
                np.max([current_max_alt, this_alt]) for current_max_alt, this_alt in zip(overall_max_alts, dist)
            ]
        if calc_moon_distance:
            # the moon doesn't move much, so just calculate the distance at the start of the night
            moon = get_body("moon", beg, observing_session.observer.location)
            moon_dist_column = f"{column_name} Min Moon Dist"
            # NOTE: for some reason it has to be moon.separation(coords) not coords.separation(moon)!!!
            dist = moon.separation(coords).value
            answer.target_list[moon_dist_column] = dist
            nightly_columns.append(moon_dist_column)
            overall_min_moon_dist = [
                np.min([current_min_dist, this_dist]) for current_min_dist, this_dist in zip(overall_min_moon_dist, dist)
            ]
    any_night_column = f"{column_prefix}Any Night"
    answer.target_list[any_night_column] = overall_any_night
    every_night_column = f"{column_prefix}Every Night"
    answer.target_list[every_night_column] = overall_every_night
    main_columns = [any_night_column, every_night_column]
    if calc_max_altitude:
        max_alt_column = f"{column_prefix}Max Alt"
        answer.target_list[max_alt_column] = overall_max_alts
        main_columns.append(max_alt_column)
    if calc_moon_distance:
        moon_dist_column = f"{column_prefix}Min Moon Dist"
        answer.target_list[moon_dist_column] = overall_min_moon_dist
        main_columns.append(moon_dist_column)
    answer.column_groups[column_prefix.strip()] = (main_columns, nightly_columns)
    moon_phases = pd.DataFrame()
    moon_phases["Time"] = [beg.iso[:10] for beg, _ in observing_session.observing_segments]
    moon_phases["Phase"] = [ap.moon_illumination(t) for t, _ in observing_session.observing_segments]
    answer.addOther("Lunar Phases", moon_phases)
    return answer


def filter_targets(
    tl: TargetList,
    criteria=True,
    inverse: bool = False,
    **kwargs,
) -> TargetList:
    answer = tl.copy()
    if inverse:
        answer.target_list = answer.target_list[~criteria(answer.target_list)]
    else:
        answer.target_list = answer.target_list[criteria(answer.target_list)]
    return answer
