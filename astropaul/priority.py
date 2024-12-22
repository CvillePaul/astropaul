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


def calculate_altitude_priority(
    pl: PriorityList, altitude_categories: CategoryTable, min_nonzero_time: u.Quantity = 0 * u.hour
) -> None:
    for _, row in pl.target_list.target_list.iterrows():
        target = row["Target Name"]
        coord = SkyCoord(ra=row["ra"], dec=row["dec"], unit=u.deg)
        for segment_table in pl.target_tables[target]:
            times = Time(segment_table.index)
            altitudes = pl.session.observer.altaz(times, coord).alt.value
            segment_table["Altitude Value"] = altitudes
            priorities = [pick_category(altitude_categories, altitude) for altitude in altitudes]
            num_nonzero = len(list(filter(lambda x: x > 0, priorities)))
            if pl.interval * num_nonzero > min_nonzero_time:
                segment_table["Altitude Priority"] = priorities
            else:
                segment_table["Altitude Priority"] = 0


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
    out_eclipse_name:str = "Not in Eclipse",
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
    else:
        raise ValueError("No target tables are present in PriorityList")

    for target_name, segment_tables in pl.target_tables.items():
        target_events = phase_events[phase_events["Target Name"] == target_name]
        target_systems = [system for system, _ in target_events.groupby("System")]
        eclipse_columns = [f"System {system} Eclipse" for system in target_systems]
        problem_columns = [f"System {system} Problems" for system in target_systems]
        for segment_table in segment_tables:
            segment_beg = Time(segment_table.index[0]).jd
            segment_end = Time(segment_table.index[-1]).jd
            # add columns for each system of this target
            for eclipse_column, problem_column in zip(eclipse_columns, problem_columns):
                segment_table[eclipse_column] = ""
                segment_table[problem_column] = ""
            # figure out all the eclipses that happen for this target during this segment
            target_eclipses = []  # list of (beg, end) jd values of eclipses.  nan for beg or end indicates partial eclipse
            for (system, member), group in target_events.groupby(["System", "Member"]):
                if group.empty or len(group) < 3:
                    continue  # skip this target if it doesn't have enough events to fulfill the full eclipse pattern
                mask = (group["Event JD"] >= segment_beg) & (group["Event JD"] <= segment_end)
                segment_events = group[mask].sort_values("Event JD")
                eclipse_beg = float("nan")
                eclipse_under_construction = False
                for _, event_name, event_jd, phase in segment_events[["Event", "Event JD", "Phase"]].itertuples():
                    if event_name == in_eclipse_name:
                        eclipse_under_construction = True
                        if phase == phase:
                            eclipse_beg = event_jd  # eclipse started during this segment
                        else:
                            eclipse_beg = float("nan")  # eclipse already in progress when segment began
                    elif event_name == out_eclipse_name and eclipse_under_construction:
                        target_eclipses.append((system, member, eclipse_beg, event_jd))
                        eclipse_under_construction = False
                if eclipse_under_construction:
                    target_eclipses.append(
                        (system, member, eclipse_beg, float("nan"))
                    )  # eclipse didn't finish before segment end
            segment_table["Eclipse Priority"] = no_eclipse_score
            if len(target_eclipses) == 0:
                continue
            # for each eclipse, calc the segment index on which the beg & end occur
            eclipse_ranges = []
            for _, _, beg, end in target_eclipses:
                beg_dt = Time(beg, format="jd").to_datetime() if beg == beg else Time(segment_beg, format="jd").to_datetime()
                end_dt = Time(end, format="jd").to_datetime() if end == end else Time(segment_end, format="jd").to_datetime()
                beg_index = segment_table.index[segment_table.index <= beg_dt].max()
                end_index = segment_table.index[segment_table.index >= end_dt].min()
                eclipse_ranges.append((beg_index, end_index))
            # check if there are problems with any eclipses
            eclipse_problems = {index: [] for index in range(len(target_eclipses))}
            # flag all partial eclipses
            for i, (_, _, beg, end) in enumerate(target_eclipses):
                if beg != beg:
                    eclipse_problems[i].append("No ingress")
                elif end != end:
                    eclipse_problems[i].append("No egress")
            # check if two eclipses overlap or one surrounds the other
            for i in range(len(target_eclipses)):
                system_1, _, beg_1, end_1 = target_eclipses[i]
                for j in range(i + 1, len(target_eclipses)):
                    system_2, _, beg_2, end_2 = target_eclipses[j]
                    if beg_1 <= beg_2 <= end_1 or beg_1 <= end_2 <= end_1 or beg_2 <= beg_1 <= end_2 or beg_2 <= end_1 <= end_2:
                        eclipse_problems[i].append(f"{system_2} Overlap")
                        eclipse_problems[j].append(f"{system_1} Overlap")
            # check if altitude ok during all of each eclipse
            for i, (beg_index, end_index) in enumerate(eclipse_ranges):
                if segment_table.loc[beg_index:end_index, altitude_column].min() < min_altitude.value:
                    eclipse_problems[i].append("Low")
            # mark appropriate rows for each eclipse
            for (system, member, _, _), (problems), (beg_index, end_index) in zip(
                target_eclipses, eclipse_problems.values(), eclipse_ranges
            ):
                segment_table.loc[beg_index:end_index, f"System {system} Eclipse"] = f"{system}{member}"
                segment_table.loc[beg_index:end_index, f"System {system} Problems"] = ", ".join(problems)
            # finally, determine priority based on absence of problems
            mask = (segment_table[eclipse_columns].ne("").any(axis=1)) & (  # require one or more eclipses for this sub-segment
                segment_table[problem_columns].eq("").all(axis=1)
            )  # require no problems for any eclipse in this sub-segment
            segment_table.loc[mask, "Eclipse Priority"] = 1.0


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
        if len(aggregate_table.columns) == 0:
            continue  # don't add tables with no columns
        segment_target_list = target_list[target_list["Target Name"].isin(aggregate_table.columns)].sort_values("ra")
        target_names = segment_target_list["Target Name"].tolist()
        # calculate LST at beginning of this observing segment
        lst = pl.session.observer.local_sidereal_time(sub_segments[0][0]).to(u.deg).value
        first_ra_of_night = (lst - 60) % 360  # lst is the RA at 90 deg alt, we want targets at 30 deg alt
        # now "pivot" the table so the first column is the RA
        i = 0
        for _, target in segment_target_list.iterrows():
            if target["ra"] > first_ra_of_night:
                break
            i += 1
        target_names = target_names[i:] + target_names[:i]  # "rotate" the list so first ra is one > lst
        pl.numerical_priorities.append(aggregate_table[target_names])
