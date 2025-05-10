import collections
from functools import partial
import inspect
from sqlite3 import Connection
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, Angle, get_body
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd

from astropaul.csv_loader.csv2sql import db_style_to_string
from astropaul.observing import ObservingSession
import astropaul.phase as ph

import __main__


class TargetList:
    def __init__(
        self,
        name: str = "Target List",
        target_list: pd.DataFrame = None,
        list_criteria: list[str] = None,
        column_groups: dict[str, tuple[list[str], list[str]]] = None,
        other_lists: dict[str, pd.DataFrame] = {},
    ):
        self.name = name
        self.target_list = target_list if target_list is not None else pd.DataFrame()
        self.list_criteria = list_criteria or []
        self.column_groups = column_groups or collections.defaultdict(list)
        self.other_lists = other_lists

    def copy(self) -> "TargetList":
        answer = TargetList(
            name=self.name,
            target_list=self.target_list.copy(),
            list_criteria=self.list_criteria.copy(),
            column_groups=self.column_groups.copy(),
            other_lists={key: value.copy() for key, value in self.other_lists.items()},
        )
        return answer

    def add_columns(self, columns: dict[str, Any] = {}):
        pass

    def add_other(self, name: str, other: pd.DataFrame):
        self.other_lists[name] = other

    def summarize(self) -> str:
        answer = f"Name: {self.name}\n"
        answer += "Criteria\n"
        if len(self.list_criteria) > 0:
            for criterion in self.list_criteria:
                answer += f"    {criterion}\n"
        else:
            answer += "    (none)\n"
        target_types = collections.Counter(self.target_list["Target Type"])
        answer += f"{len(self.target_list)} targets:\n"
        for type, count in target_types.items():
            answer += f"    {count:4d} {type}\n"
        answer += "Column Count (primary, secondary):\n"
        for name, (primary, secondary) in self.column_groups.items():
            answer += f"    {name}: ({len(primary)}, {len(secondary)})\n"
        answer += "Associated tables:\n"
        for name, other in self.other_lists.items():
            answer += f"    {len(other):4d} rows, {len(other.columns):2d} columns: {name}\n"
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
        answer = first.copy()
        if answer.target_list.empty or answer.target_list is None:
            answer.target_list = target_list
        else:
            answer.target_list = answer.target_list.join(target_list)
        for name, (primary, secondary) in column_groups.items():
            if entry := answer.column_groups.get(name, None):
                answer.column_groups[name] = (entry[0] + primary, entry[1] + secondary)
            else:
                answer.column_groups[name] = (primary, secondary)
        for name, other_list in other_lists.items():
            answer.other_lists[name] = other_list
        return answer


class TargetListCreator:
    def __init__(self, name: str = "Standard", connection: Connection = None, steps: list = None, **kwargs):
        self.name = name
        self.connection = connection
        self.kwargs = kwargs.copy()
        self.steps = steps

    def calculate(self, initial_list: TargetList = None, steps=None, verbose: bool = False, **kwargs) -> TargetList:
        """
        Create a new target list by running through self.steps and returning the result
        """
        if not initial_list:
            initial_list = TargetList(name=self.name)
        else:
            initial_list = initial_list.copy()
            initial_list.name = self.name
        intermediate_tl = initial_list if initial_list else TargetList(name=self.name)
        merged_kwargs = {"connection": self.connection, **self.kwargs, **kwargs}
        if steps is None:
            steps = self.steps
        for step in steps:
            intermediate_tl = step(intermediate_tl, **merged_kwargs)
            if verbose:
                print(f"{len(intermediate_tl.target_list):4d} targets, {step=}")
        if verbose:
            print(intermediate_tl.summarize())
        return intermediate_tl


# some convenience functions
def convert_columns_to_human_style(df: pd.DataFrame) -> None:
    if index_name := df.index.name:
        df.index.name = db_style_to_string(index_name)
    df.columns = map(db_style_to_string, df)


def verify_step_requirements(tl: TargetList, required_tables: set[str] = None, required_columns: set[str] = None):
    if tl.target_list.empty:
        raise ValueError("Target List cannot be empty")
    if required_tables:
        existing_tables = set(tl.other_lists.keys())
        if not required_tables.issubset(existing_tables):
            raise ValueError(f"Required table(s) missing: {required_tables - existing_tables}")
    if required_columns:
        existing_columns = set(tl.target_list.columns)
        if not required_columns.issubset(existing_columns):
            raise ValueError(f"Required column(s) missing: {required_columns - existing_columns}")


# following are methods to be used as steps of a TargetListCreator


def add_targets(tl: TargetList = TargetList(), **kwargs) -> TargetList:
    """Add rows for each target in the database"""
    targets = pd.read_sql(
        "select * from targets;",
        kwargs["connection"],
    )
    convert_columns_to_human_style(targets)
    column_groups = {"Target": (["Target Name"], [column for column in targets.columns if column != "Target Name"])}
    return TargetList.merge(tl, targets, column_groups, {})


def concat_dataframe(tl: TargetList, other_df: pd.DataFrame, **kwargs) -> TargetList:
    return TargetList(
        name=tl.name,
        target_list=pd.concat([tl.target_list, other_df]),
        column_groups=tl.column_groups,
        other_lists=tl.other_lists,
    )


def add_lists(tl: TargetList, **kwargs) -> TargetList:
    verify_step_requirements(tl)
    conn = kwargs["connection"]
    target_lists = [name[0] for name in conn.execute("select target_list from target_lists;").fetchall()]
    answer = tl.copy()
    for target_list in target_lists:
        list_members = [
            result[0]
            for result in conn.execute(
                """
                select target_name 
                from target_list_members
                where target_list = ?
                ;""",
                [target_list],
            ).fetchall()
        ]
        column_name = f"List {target_list}"
        answer.target_list[column_name] = answer.target_list["Target Name"].isin(list_members)
    answer.column_groups["List"] = ([], [f"List {name}" for name in target_lists])
    return answer


def add_phase_events(
    tl: TargetList, observing_session: ObservingSession, phase_event_defs: list[ph.PhaseEventDef], **kwargs
) -> TargetList:
    """Calculate what phase events are in effect for each observing segment"""
    verify_step_requirements(tl, {"Ephemerides"})
    ephem_table = tl.other_lists["Ephemerides"]
    answer = tl.copy()
    ephems = {}
    phase_events = pd.DataFrame(
        columns=["Target Name", "System", "Member", "Orbit Num", "Beg JD", "Phase", "Event", "End JD", "Beg UTC", "End UTC"]
    )
    for target_name, ephem_row in ephem_table.sort_values(["System", "Member"]).iterrows():
        ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
        ephems[target_name] = ephems.get(target_name, []) + [ephem]
    for beg, end in observing_session.observing_segments:
        for target_name, ephem_list in ephems.items():
            for ephem in ephem_list:
                event_list = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg.jd, end.jd).events
                event_list[0].jd = beg.jd  # treat first event as if it had started at segment beginning
                # set end time equal to begin time of next event, or if last event to segment end
                num_events = len(event_list)
                end_jds = [0.0] * num_events
                end_utcs = [""] * num_events
                for i in range(num_events):
                    if i == num_events - 1:
                        end_jd = end.jd
                    else:
                        end_jd = event_list[i + 1].jd
                    end_jds[i] = end_jd
                    end_utcs[i] = Time(end_jd, format="jd").iso[:19]

                # populate table with all values except end times
                for event, end_jd, end_utc in zip(event_list, end_jds, end_utcs):
                    phase_events.loc[len(phase_events)] = (
                        [target_name] + event.to_list() + [end_jd, Time(event.jd, format="jd").iso[:19], end_utc]
                    )
    answer.other_lists["Phase Events"] = phase_events
    return answer


def ancillary_data_from_tess(tl: TargetList, **kwargs) -> TargetList:
    """Add fundamental things like magnitude and effective temperature if available from TESS catalog"""
    answer = tl.copy()
    tess_data = pd.read_sql(
        """select 'TIC ' || id 'Target Name', pmra 'PM RA', pmdec 'PM Dec', vmag Vmag, teff Teff, plx Parallax, d Distance
        from tess_ticv8 tt;""",
        kwargs["connection"],
    )
    answer.target_list = answer.target_list.merge(tess_data, on="Target Name", how="left")
    answer.column_groups = {
        **answer.column_groups,
        "TESS Data": (["Vmag", "Teff"], ["PM RA", "PM Dec", "Distance", "Parallax"]),
    }
    return answer


def hide_cols(tl: TargetList, prefix: str, **kwargs) -> TargetList:
    answer = tl.copy()
    cols_to_remove = [col for col in tl.target_list.columns if col.startswith(prefix)]
    if len(cols_to_remove) > 0:
        answer.target_list.drop(labels=cols_to_remove, axis=1, inplace=True)
    del answer.column_groups[prefix]
    return answer


def add_observability(
    tl: TargetList,
    observing_session: ObservingSession,
    observability_threshold: tuple[u.Quantity, u.Quantity] = (30 * u.deg, 80 * u.deg),
    constraints: list[ap.Constraint] = None,
    column_prefix: str = "Observable ",
    calc_moon_distance: bool = False,
    **kwargs,
) -> TargetList:
    lo_alt, hi_alt = observability_threshold[0].to(u.deg).value, observability_threshold[1].to(u.deg).value
    if not constraints:
        constraints = [
            ap.AltitudeConstraint(lo_alt, hi_alt, True),
        ]
    coords = SkyCoord(ra=tl.target_list["RA"].values * u.deg, dec=tl.target_list["Dec"].values * u.deg)

    # calculate observability for each time segment of the observing session
    answer = tl.copy()
    answer.list_criteria.append(f"Observability threshold: {observability_threshold}")
    for constraint in constraints:
        answer.list_criteria.append(f"{constraint.__class__.__name__}: {vars(constraint)}")
    overall_any_night = np.array([False] * len(answer.target_list))
    overall_every_night = np.array([True] * len(answer.target_list))
    overall_max_alts = np.array([-90.0] * len(answer.target_list))
    overall_min_hours_observable = np.array([24] * len(answer.target_list))
    overall_min_moon_dist = np.array([180.0] * len(answer.target_list))
    nightly_columns = []
    time_resolution = 15 * u.min
    for beg, end in observing_session.observing_segments:
        time_grid = ap.utils.time_grid_from_range((beg, end), time_resolution=time_resolution)
        segment_alts = [observing_session.observer.altaz(time_grid, coord).alt.value for coord in coords]
        segment_good_alts = [
            len(target_alts[(lo_alt <= target_alts) & (target_alts <= hi_alt)]) for target_alts in segment_alts
        ]
        segment_hours_observable = (segment_good_alts * time_resolution).to(u.hour).value
        segment_max_alts = [np.max(target_alts) for target_alts in segment_alts]
        segment_observability = list(map(lambda x: x > 0, segment_good_alts))
        column_name = f"{column_prefix}{beg.iso[:10]}"
        answer.target_list[column_name] = segment_observability
        overall_any_night |= segment_observability
        overall_every_night &= segment_observability
        nightly_columns.append(column_name)
        hours_observable_column = f"{column_name} Hours Observable"
        overall_min_hours_observable = [
            np.min([current_min_hours, this_min_hours])
            for current_min_hours, this_min_hours in zip(overall_min_hours_observable, segment_hours_observable)
        ]
        answer.target_list[hours_observable_column] = segment_hours_observable
        nightly_columns.append(hours_observable_column)
        max_alt_column = f"{column_name} Max Alt"
        answer.target_list[max_alt_column] = segment_max_alts
        nightly_columns.append(max_alt_column)
        overall_max_alts = [
            np.max([current_max_alt, this_alt]) for current_max_alt, this_alt in zip(overall_max_alts, segment_max_alts)
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
    max_alt_column = f"{column_prefix}Max Alt"
    answer.target_list[max_alt_column] = overall_max_alts
    main_columns.append(max_alt_column)
    min_hours_observable_column = f"{column_prefix}Hours Min"
    answer.target_list[min_hours_observable_column] = overall_min_hours_observable
    main_columns.append(min_hours_observable_column)
    if calc_moon_distance:
        moon_dist_column = f"{column_prefix}Min Moon Dist"
        answer.target_list[moon_dist_column] = overall_min_moon_dist
        main_columns.append(moon_dist_column)
    answer.column_groups[column_prefix.strip()] = (main_columns, nightly_columns)
    moon_phases = pd.DataFrame()
    moon_phases["Time"] = [beg.iso[:10] for beg, _ in observing_session.observing_segments]
    moon_phases["Phase"] = [ap.moon_illumination(t) for t, _ in observing_session.observing_segments]
    answer.add_other("Lunar Phases", moon_phases)
    return answer


def filter_targets(
    tl: TargetList,
    criteria,
    inverse: bool = False,
    **kwargs,
) -> TargetList:
    answer = tl.copy()
    code = inspect.getsource(criteria).strip()
    if (prefix := code.find("criteria=")) > 0:
        code = code[prefix + 9 :]  # try to isolate the code of the lambda function passed in
    for variable in criteria.__code__.co_names:
        if value := __main__.__dict__.get(variable, None):
            code = code.replace(variable, repr(value))
    if inverse:
        code = "Inverse of: " + code
        answer.target_list = answer.target_list[~criteria(answer.target_list)]
    else:
        answer.target_list = answer.target_list[criteria(answer.target_list)]
    answer.list_criteria.append(code)
    return answer


def add_dssi_phase(tl: TargetList, phase_event_defs: list[ph.PhaseEventDef], **kwargs) -> TargetList:
    """Calculate what PhaseEventDef was most in effect during a DSSI observation"""
    verify_step_requirements(tl, {"DSSI Observations", "Ephemerides"})
    answer = tl.copy()
    ephem_table = tl.other_lists["Ephemerides"]
    dssi_phases = pd.DataFrame(columns=["Target Name", "DSSI Session", "JD Mid", "UTC Mid", "System", "Member", "State"])
    for _, dssi in tl.other_lists["DSSI Observations"].iterrows():
        target_name = dssi["Target Name"]
        ephem_rows = ephem_table[ephem_table["Target Name"] == target_name]
        if ephem_rows.empty:
            continue
        beg, end = dssi["Start JD"], dssi["End JD"]
        for _, ephem_row in ephem_rows.sort_values(["System", "Member"]).iterrows():
            ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
            state = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg, end).calc_longest_span(beg, end)
            dssi_phases.loc[len(dssi_phases)] = [
                target_name,
                int(dssi["DSSI Session"]),
                dssi["Mid JD"],
                dssi["Mid UTC"],
                ephem_row["System"],
                ephem_row["Member"],
                state,
            ]
    if not dssi_phases.empty:
        answer.other_lists["DSSI Phases"] = dssi_phases
    return answer


def add_system_configuration(
    tl: TargetList, table_name: str, time_column: str, time_format: str = "jd", **kwargs
) -> TargetList:
    """For the specified table, add a pair of columns indicating the state of each system of the target.
    For these calculations, the specified column gives the time, and must be in the specified format.
    For a system in eclipse, show what member and the percentage of the time between ingress and egress.
    For a system not in eclipse, show the percent of phase, where 0.0 is mid eclipse of the a member."""
    verify_step_requirements(tl, {table_name})
    table = tl.other_lists[table_name]
    if not time_column in table.columns:
        raise ValueError(f"Table {table_name} does not have a column {time_column}")
    answer = tl.copy()
    all_targets = set(table["Target Name"].unique())
    # determine max number of systems we need to handle, max duration
    ephem = tl.other_lists["Ephemerides"]
    ephem = ephem[ephem["Target Name"].isin(all_targets)]
    systems = sorted(ephem["System"].unique())
    # add the desired columns to the table
    phase_event_defs = [
        ph.PhaseEventDef("Not in Eclipse", partial(ph.calc_time_of_gress, ingress=False)),
        ph.PhaseEventDef("Eclipse", partial(ph.calc_time_of_gress, ingress=True)),
    ]
    system_eclipses, system_durations = {}, {}
    for system in systems:
        system_eclipses[system] = []
        system_durations[system] = []
    for _, row in table.iterrows():
        target_name = row["Target Name"]
        observation_time = Time(row[time_column], format=time_format).jd
        target_ephem = ephem[ephem["Target Name"] == target_name]
        for system in systems:
            system_ephem = target_ephem[target_ephem["System"] == system]
            if system_ephem.empty:
                system_eclipses[system].append("")
                system_durations[system].append(float("nan"))
                continue
            eclipse_found = False
            for _, ephem_row in system_ephem.sort_values("Member").iterrows():
                member_ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
                if member_ephem.duration != member_ephem.duration:
                    continue
                event_list = ph.PhaseEventList.calc_phase_events(
                    member_ephem, phase_event_defs, observation_time, observation_time + member_ephem.period
                )
                if len(event_list.events) < 2:
                    raise ValueError(f"Unexpected event list {event_list}")
                if event_list.events[0].type == "Eclipse":
                    system_eclipses[system].append(f"{system}{ephem_row["Member"]}")
                    duration_percent = 1 - (event_list.events[1].jd - observation_time) / member_ephem.duration
                    system_durations[system].append(duration_percent)
                    eclipse_found = True
                    break
            if not eclipse_found:
                system_eclipses[system].append("")
                system_durations[system].append(float("nan"))
    for system in systems:
        table[f"System {system} Eclipse"] = system_eclipses[system]
        table[f"System {system} Duration Percent"] = system_durations[system]
    answer.other_lists[table_name] = table
    return answer


def add_side_status(
    tl: TargetList, phase_event_defs: list[ph.PhaseEventDef], side_state: str = "Eclipse", **kwargs
) -> TargetList:
    """Examine all speckle observations and determine which ones occurred during eclipse.
    `side_state` specifies the name returned by one or more `phase_event_defs` that indicate an eclipse."""
    verify_step_requirements(tl, {"Ephemerides", "DSSI Observations"})
    answer = tl.copy()
    ephem_table = tl.other_lists["Ephemerides"]
    side_observations = pd.DataFrame(
        columns=["Target Name", "DSSI Session", "Mid JD", "Mid UTC", "System", "Member", "SIDE Type"]
    )
    for target_name, dssi in tl.other_lists["DSSI Observations"].iterrows():
        ephem_rows = ephem_table[ephem_table["Target Name"] == target_name]
        if ephem_rows.empty:
            continue
        beg, end = dssi["Start JD"], dssi["End JD"]
        for _, ephem_row in ephem_rows.sort_values(["System", "Member"]).iterrows():
            ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
            state = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg, end).calc_longest_span(beg, end)
            if state == side_state:
                if ephem.duration == ephem.duration:
                    side_type = "Exact"
                else:
                    side_type = "Synthetic Duration"
                side_observations.loc[len(side_observations)] = [
                    target_name,
                    int(dssi["DSSI Session"]),
                    dssi["Mid JD"],
                    dssi["Mid UTC"],
                    ephem_row["System"],
                    ephem_row["Member"],
                    side_type,
                    side_type,
                ]
    if not side_observations.empty:
        answer.other_lists["SIDE Observations"] = side_observations.sort_values(["Target Name", "System", "Member", "JD Mid"])
        side_counts = side_observations.groupby("Target Name").size()
        side_counts.name = "Num SIDE"
        answer.target_list = answer.target_list.merge(side_counts, on="Target Name", how="left")
        answer.target_list["Num SIDE"] = answer.target_list["Num SIDE"].fillna(0).astype(int)
    return answer


def add_rv_status(tl: TargetList, phase_event_defs: list[ph.PhaseEventDef], **kwargs) -> TargetList:
    """Calculate what PhaseEventDef was most in effect during a PEPSI spectrum.
    The status of each system is determined by the status of its first component with non-nan ephemeris values.
    These system statuses are then concatenated, separated by '|' characters, to get the final status for the target.
    """
    verify_step_requirements(tl, {"PEPSI Observations", "Ephemerides"})
    answer = tl.copy()
    ephem_table = tl.other_lists["Ephemerides"]
    pepsi_table = tl.other_lists["PEPSI Observations"]
    pepsi_phases = pd.DataFrame(columns=["Target Name", "Start JD", "Start UTC", "State"])
    for _, pepsi in pepsi_table.iterrows():
        ephem_rows = ephem_table[ephem_table["Target Name"] == pepsi["Target Name"]]
        if ephem_rows.empty:
            continue
        exposure_days = pepsi["Exposure Time"] / 24 / 3600  # exposure time is in seconds instead of days
        beg = pepsi["Mid JD"] - exposure_days / 2
        end = beg + exposure_days
        system_states = []
        for _, system_ephems in ephem_rows.groupby("System"):
            system_state = None
            for _, component_ephem in system_ephems.iterrows():
                ephem = ph.Ephemeris.from_dataframe_row(component_ephem)
                if ephem.t0 == ephem.t0 and ephem.period == ephem.period and ephem.duration == ephem.duration:
                    system_state = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg, end).calc_longest_span(
                        beg, end
                    )
                    break
            if system_state:
                system_states.append(system_state)
            else:
                system_states.append("?")

        pepsi_phases.loc[len(pepsi_phases)] = [
            pepsi["Target Name"],
            beg,
            Time(beg, format="jd").iso[:18],
            "|".join(system_states),
        ]
    if not pepsi_phases.empty:
        answer.other_lists["PEPSI RV Status"] = pepsi_phases
    return answer


def add_tess_sectors(tl: TargetList, **kwargs) -> TargetList:
    from astroquery.mast import Tesscut

    coords = SkyCoord(ra=tl.target_list["ra"], dec=tl.target_list["dec"], unit=u.deg)
    answer = tl.copy()
    answer.target_list["TESS Sectors"] = [
        ", ".join(Tesscut.get_sectors(coordinates=coord)["sector"].astype(str)) for coord in coords
    ]
    answer.target_list["TESS Sectors"] = [
        ", ".join(Tesscut.get_sectors(coordinates=coord)["sector"].astype(str)) for coord in coords
    ]


def add_catalogs(tl: TargetList, **kwargs) -> TargetList:
    answer = tl.copy()
    conn = kwargs["connection"]
    all_targets = list(tl.target_list["Target Name"].unique())
    catalog_members = pd.read_sql("select target_name, catalog_name, catalog_id from catalog_members;", conn)
    tess_catalog_ids = pd.read_sql(
        """
        select cm.target_name, tt.gaia 'GAIA DR2', hip, twomass
        from catalog_members cm
        join tess_ticv8 tt on tt.id = cm.catalog_id
        where cm.catalog_name = 'TESS TICv8';
   
    """,
        conn,
    )
    catalog_members_from_tess = tess_catalog_ids.melt(id_vars="target_name", var_name="catalog_name", value_name="catalog_id")
    catalog_members = pd.concat([catalog_members, catalog_members_from_tess], ignore_index=True)
    convert_columns_to_human_style(catalog_members)
    catalog_members = catalog_members[catalog_members["Target Name"].isin(all_targets)]
    answer.other_lists["Catalog Membership"] = catalog_members
    return answer


def add_database_table(
    tl: TargetList, table_name: str, add_count: bool = True, **kwargs
) -> TargetList:
    verify_step_requirements(tl)
    answer = tl.copy()
    table_contents = pd.read_sql(f"select * from {table_name};", kwargs["connection"])
    convert_columns_to_human_style(table_contents)
    human_name = db_style_to_string(table_name)
    if "Target Name" in table_contents.columns:
        all_targets = list(tl.target_list["Target Name"].unique())
        table_contents = table_contents[table_contents["Target Name"].isin(all_targets)]
        if add_count:
            count_column = f"Num {human_name}"
            counts = table_contents["Target Name"].value_counts()
            answer.target_list[count_column] = answer.target_list["Target Name"].map(counts).fillna(0).astype(int)
            primary_columns, secondary_columns = answer.column_groups.get("Count", ([], []))
            primary_columns.append(count_column)
            answer.column_groups["Count"] = (primary_columns, secondary_columns)
    answer.other_lists[human_name] = table_contents
    return answer
