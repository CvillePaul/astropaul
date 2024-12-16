import operator
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, get_body
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd

import astropaul.targetlistcreator as tlc
import astropaul.phase as ph


class PriorityList:
    def __init__(
        self,
        target_list: tlc.TargetList,
        session: tlc.ObservingSession,
        interval: u.Quantity = 1 * u.hour,
    ):
        self.target_list = target_list
        self.session = session
        self.interval = interval
        self.targets = target_list.target_list["Target Name"]  # TODO: remove hard coded column name
        self.segments = session.calc_subsegments(interval)
        self.target_tables = {}  # key=target name, value = list[dataframe], one per observing segment, indexed by subsegment
        for target in self.targets:
            tables = []
            for segment in self.segments:
                table = pd.DataFrame(index=[beg.to_datetime() for beg, _ in segment])
                tables.append(table)
            self.target_tables[target] = tables
        self.numerical_priorities: list[pd.DataFrame] = None
        self.categorical_priorities: list[pd.DataFrame] = None
        self.category_bins = None
        self.category_labels = None

    def categorize_priorities(self, bins: list[float] = None, labels: list[str] = None) -> None:
        if not bins or not labels:
            raise ValueError("No argument can be None")
        if len(bins) != len(labels) + 1:
            raise ValueError("List of bins must be one element larger than list of labels")
        if not self.numerical_priorities:
            raise ValueError("Overall numerical priorities have not been calculated yet")
        self.category_bins = bins
        self.category_labels = labels
        self.categorical_priorities = []
        for numerical_priority in self.numerical_priorities:
            category_table = numerical_priority.apply(
                lambda col: pd.cut(col, bins=bins, labels=labels, ordered=False, include_lowest=True)
            )
            self.categorical_priorities.append(category_table)


CategoryTable = list[tuple[tuple[Any, Any], Any]]


def pick_category(categories: CategoryTable, value: float) -> Any:
    for (beg, end), result in categories:
        if beg <= value < end:
            return result


def calculate_moon_priority(
    pl: PriorityList, illumination_categories: CategoryTable = None, dist_categories: CategoryTable = None
) -> None:
    for _, row in pl.target_list.target_list.iterrows():
        target = row["Target Name"]
        coord = SkyCoord(ra=row["ra"], dec=row["dec"], unit=u.deg)
        for segment_table in pl.target_tables[target]:
            times = Time(segment_table.index)
            illuminations = ap.moon_illumination(times)
            illumination_names = [pick_category(illumination_categories, illumination) for illumination in illuminations]
            moon = get_body("moon", times, pl.session.observer.location)  # .transform_to(GCRS(obstime=times))
            distances = moon.separation(coord).value
            priorities = [
                pick_category(dist_categories[illumination_name], distance)
                for illumination_name, distance in zip(illumination_names, distances)
            ]
            segment_table["Moon Illumination"] = illuminations
            segment_table["Moon Illumination Name"] = illumination_names
            segment_table["Moon Distance"] = distances
            segment_table["Moon Priority"] = priorities


def calculate_altitude_priority(pl: PriorityList, altitude_categories: CategoryTable) -> None:
    for _, row in pl.target_list.target_list.iterrows():
        target = row["Target Name"]
        coord = SkyCoord(ra=row["ra"], dec=row["dec"], unit=u.deg)
        for segment_table in pl.target_tables[target]:
            times = Time(segment_table.index)
            altitudes = pl.session.observer.altaz(times, coord).alt.value
            priorities = [pick_category(altitude_categories, altitude) for altitude in altitudes]
            segment_table["Altitude Value"] = altitudes
            segment_table["Altitude Priority"] = priorities


def calculate_list_priority(pl: PriorityList, list_name: str, invert: bool = False, false_value: float = 0.2) -> None:
    for _, row in pl.target_list.target_list.iterrows():
        target = row["Target Name"]
        list_value = row[f"List {list_name}"]
        if list_value != list_value:  # set NaN values to False
            list_value = False
        priority = 1 if (list_value ^ invert) else false_value
        for table in pl.target_tables[target]:
            table[f"{list_name} Priority"] = priority


def calculate_prior_observation_priority(pl: PriorityList, prior_observation_categories: CategoryTable) -> None:
    pass


def calculate_eclipse_priority(
    pl: PriorityList,
    in_eclipse_name: str = "Eclipse",
    out_eclipse_name="Not in Eclipse",
    no_eclipse_score: float = 0.2,
    min_altitude: u.Quantity = 30 * u.deg,
) -> None:
    """Favor targets that undergo a full ingress/eclipse/egress during an observing segment"""
    required_table = "Phase Events"
    phase_events = pl.target_list.other_lists.get(required_table, pd.DataFrame())
    if phase_events.empty:
        raise ValueError(f"TargetList does not have a {required_table} table")
    altitude_column = "Altitude Value"
    if pl.target_tables:
        if not altitude_column in next(iter(pl.target_tables.values()))[0].columns:
            raise ValueError(f"Required column {altitude_column} not in existing target tables")
    eclipse_pattern = [out_eclipse_name, in_eclipse_name, out_eclipse_name]
    for _, row in pl.target_list.target_list.iterrows():
        target_name = row["Target Name"]
        # figure out all the full eclipses that happen for this target
        target_events = phase_events[phase_events["Target Name"] == target_name]
        target_eclipses = []
        for (system, member), group in target_events.groupby(["System", "Member"]):
            if group.empty or len(group) < len(eclipse_pattern):
                continue
            times = []
            for _, event in group.sort_values("Event JD").iterrows():
                if event["Event"] == eclipse_pattern[len(times)]:
                    times.append(event["Event JD"])
                    if len(times) == len(eclipse_pattern):
                        break
                    continue
            if len(times) == len(eclipse_pattern):
                beg = Time(times[1], format="jd").to_datetime()
                end = Time(times[2], format="jd").to_datetime()
                target_eclipses.append((system, member, beg, end))
        # if two eclipses overlap or one surrounds the other, remove both
        for i in range(len(target_eclipses)):
            system_1, _, beg_1, end_1 = target_eclipses[i]
            if system_1 == "":
                continue
            for j in range(i + 1, len(target_eclipses)):
                system_2, _, beg_2, end_2 = target_eclipses[j]
                if system_2 == "":
                    continue
                if beg_1 <= beg_2 <= end_1 or beg_1 <= end_2 <= end_1:
                    target_eclipses[i] = ("", "", beg_1, end_1)
                    target_eclipses[j] = ("", "", beg_2, end_2)
        target_eclipses = [x for x in target_eclipses if x[0] != ""]
        # add columns to the segment table indicating the member in eclipse for each system
        system_columns = {system: f"Full {system} Eclipse" for system, _ in target_events.groupby("System")}
        if len(system_columns) == 0:
            continue  # skip if none of the systems ever have an eclipse
        for segment_table in pl.target_tables[target_name]:
            for col_name in system_columns.values():
                segment_table[col_name] = ""
            segment_beg = segment_table.index[0]
            segment_end = segment_table.index[-1]
            for system, member, beg, end in target_eclipses:
                if segment_beg <= beg and end <= segment_end:
                    beg_index = segment_table.index[segment_table.index <= beg].max()
                    end_index = segment_table.index[segment_table.index >= end].min()
                    if segment_table.loc[beg_index:end_index, altitude_column].min() >= min_altitude.value:
                        segment_table.loc[beg_index:end_index, system_columns[system]] = f"{system}{member}"
            # set the eclipse priority based on the eclipse status of each system's column
            segment_table["Eclipse Priority"] = no_eclipse_score
            segment_table.loc[segment_table[system_columns.values()].sum(axis=1).str.len() == 2, "Eclipse Priority"] = 1.0


def calculate_phase_priority(pl: PriorityList, phase_defs: list[ph.PhaseEventDef], phase_categories: dict[str, float]) -> None:
    ephem_table = pl.target_list.other_lists["Ephem"]
    min_priority = min(phase_categories, key=operator.itemgetter(1))
    for target, table_list in pl.target_tables.items():
        ephem_rows = ephem_table.loc[target]
        ephem_rows = ephem_rows[ephem_rows["Ephem Member"] == "a"]
        if ephem_rows.empty or len(ephem_rows) < 2:
            # flunk anything that doesn't have at least two binaries
            for table in table_list:
                table["Phase Priority"] = min_priority
            continue
        for _, row in ephem_rows.iterrows():
            ephem = ph.Ephemeris.from_dataframe_row(row)
            session_beg, session_end = pl.session.time_range
            event_list = ph.PhaseEventList.calc_phase_events(ephem, phase_defs, session_beg.jd, session_end.jd)
            # make columns for each system for each observing segment
            for table, segments in zip(table_list, pl.segments):
                states = [event_list.calc_longest_span(segment_beg.jd, segment_end.jd) for segment_beg, segment_end in segments]
                table[f"Phase System {ephem.system} State"] = states
        # make the overall priority column that considers all available systems
        system_cols = [f"Phase System {system} State" for system in ephem_rows["Ephem System"]]
        for table in table_list:
            overall_priorities = []
            for _, row in table.iterrows():
                whole_state = "|".join(sorted(row[system_cols]))
                overall_priorities.append(phase_categories[whole_state])
            table["Phase Priority"] = overall_priorities


def calculate_overall_priority(pl: PriorityList) -> None:
    for table_list in pl.target_tables.values():
        for table in table_list:
            priority_columns = [col for col in table.columns if col.endswith(" Priority")]
            table["Overall Priority"] = 1
            for col in priority_columns:
                table["Overall Priority"] *= table[col]


def aggregate_target_priorities(pl: PriorityList, skip_column_threshold: float = 0) -> None:
    aggregate_tables = []
    for segment in pl.segments:
        aggregate_tables.append(pd.DataFrame(index=[beg.to_datetime() for beg, _ in segment]))
    # add a column to each table for each target's overall priority
    for target, target_tables in pl.target_tables.items():
        for target_table, aggregate_table in zip(target_tables, aggregate_tables):
            if np.max(target_table["Overall Priority"]) >= skip_column_threshold:
                aggregate_table[target] = target_table["Overall Priority"]
    pl.numerical_priorities = []
    target_list = pl.target_list.target_list
    # each observing segment might have a different list of targets, so sorting by RA needs to happen per observing segment
    for sub_segments, aggregate_table in zip(pl.segments, aggregate_tables):
        segment_target_list = target_list[target_list["Target Name"].isin(aggregate_table.columns)].sort_values("ra")
        target_names = segment_target_list["Target Name"].tolist()
        # calculate LST at beginning of this observing segment
        lst = pl.session.observer.local_sidereal_time(sub_segments[0][0]).to(u.deg).value
        first_ra_of_night = lst - 60  # lst is the RA at 90 deg alt, we want targets at 30 deg alt
        # now "pivot" the table so the first column is the RA
        i = 0
        for _, target in segment_target_list.iterrows():
            if target["ra"] > first_ra_of_night:
                break
            i += 1
        target_names = target_names[i:] + target_names[:i]  # "rotate" the list so first ra is one > lst
        pl.numerical_priorities.append(aggregate_table[target_names])
