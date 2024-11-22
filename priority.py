from datetime import datetime
import operator
import platform
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, get_body, GCRS
from astropy.time import Time
import astropy.units as u
import itables
import numpy as np
import pandas as pd

import astropaul.targetlistcreator as tlc
import astropaul.phase as ph


def dataframe_to_datatable(table: pd.DataFrame, table_name: str = "table", caption: str = "", table_options: dict = {}):
    default_options = {
        "connected": True,
        "paging": False,
        "maxBytes": 0,
        "maxColumns": 0,
        "autoWidth": False,
        "layout": {"topStart": None, "topEnd": None},
        "classes": "display nowrap compact cell-border",
    }
    if caption == "":
        caption = f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}"
    html = itables.to_html_datatable(df=table, caption=caption, table_id=table_name, **{**default_options, **table_options})
    html += f"""
        <style>
        #{table_name} th {{
            white-space: normal;
            word-wrap: break-word;
            text-align: center;
        }}
        caption {{
            text-align: left;
            font-style: italic;
        }}
        </style>
        """
    return html


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
        self.categorized_priorities: list[pd.DataFrame] = None

    def categorize_priorities(self, bins: list[float] = None, labels: list[str] = None) -> None:
        if not bins or not labels:
            raise ValueError("No argument can be None")
        if len(bins) != len(labels) + 1:
            raise ValueError("List of bins must be one element larger than list of labels")
        if not self.numerical_priorities:
            raise ValueError("Overall numerical priorities have not been calculated yet")
        self.categorized_priorities = []
        for numerical_priority in self.numerical_priorities:
            self.categorized_priorities.append(
                numerical_priority.apply(lambda col: pd.cut(col, bins=bins, labels=labels, ordered=False, include_lowest=True))
            )

    def to_html(self, file_base: str = "Priorities"):
        for target, target_table_list in self.target_tables.items():
            for i, target_table in enumerate(target_table_list):
                file = (f"Target {file_base} {target} {i:02d}.html"
                    # f"{target_table.index[0].strftime("%Y-%m-%d %H_%M")} to "
                    # f"{target_table.index[-1].strftime("%Y-%m-%d %H_%M")}.html"
                )
                with open(file, "w") as f:
                    f.write(dataframe_to_datatable(target_table, "targetPriority"))
        if self.numerical_priorities:
            for i, priority_table in enumerate(self.numerical_priorities):
                numerical_html = dataframe_to_datatable(priority_table, "numericalPriority")
                with open(f"Numerical {file_base} {i:02d}.html", "w") as f:
                    f.write(numerical_html)
        if self.categorized_priorities:
            for i, categories_table in enumerate(self.categorized_priorities):
                categorized_html = dataframe_to_datatable(categories_table, "categorizedPriority")
                with open(f"Categorized {file_base} {i:02d}.html", "w") as f:
                    f.write(categorized_html)


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


def calculate_phase_priority(pl: PriorityList, phase_defs: list[ph.PhaseEventDef], phase_categories: dict[str, float]) -> None:
    ephem_table = pl.target_list.other_lists["Ephem"]
    min_priority = min(phase_categories, key=operator.itemgetter(1))
    for target, table_list in pl.target_tables.items():
        ephem_rows = ephem_table[(ephem_table["Ephem Name"] == target) & (ephem_table["Ephem Member"] == "a")]
        if len(ephem_rows) < 2:
            # flunk anything that doesn't have at least two binaries
            for table in table_list:
                table["Phase Priority"] = min_priority
            continue
        for _, row in ephem_rows.iterrows():
            ephem = ph.Ephemeris(*row[["Ephem System", "Ephem Member", "Ephem T0", "Ephem Period", "Ephem Duration"]])
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


def calculate_overall_priority(pl: PriorityList, weightings: dict[str, float]) -> None:
    # uncomment lines in order to switch to weighted sum priorities instead of multiplicative priorities
    # norm_factor = np.sum([weight for weight in weightings.values()])
    for table_list in pl.target_tables.values():
        for table in table_list:
            table["Overall Priority"] = 1  # 0
            for col, weight in weightings.items():
                table["Overall Priority"] *= table[f"{col} Priority"] * weight
                # table["Overall Priority"] += table[f"{col} Priority"]  * weight
            # table["Overall Priority"] /= norm_factor


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
        this_target_list = target_list[target_list["Target Name"].isin(aggregate_table.columns)].sort_values("ra")
        target_names = this_target_list["Target Name"].tolist()
        # calculate LST at beginning of this observing segment
        lst = pl.session.observer.local_sidereal_time(sub_segments[0][0]).to(u.deg).value
        # now "pivot" the table so the first column is the RA 
        i = next(i for i, ra in enumerate(this_target_list["ra"].tolist()) if ra > lst) # find first element with ra > lst
        target_names = target_names[i:] + target_names[:i] # "rotate" the list so first ra is one > lst
        pl.numerical_priorities.append(aggregate_table[target_names])
