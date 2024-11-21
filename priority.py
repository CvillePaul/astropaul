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
        self.tables = {}  # key=target name, value = list[dataframe], one per observing segment, indexed by observing subsegment
        for target in self.targets:
            tables = []
            for segment in self.segments:
                table = pd.DataFrame(index=[beg.to_datetime() for beg, _ in segment])
                tables.append(table)
            self.tables[target] = tables
        self.target_priorities: list[pd.DataFrame] = None


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
        for segment_table in pl.tables[target]:
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
        for segment_table in pl.tables[target]:
            times = Time(segment_table.index)
            altitudes = pl.session.observer.altaz(times, coord).alt.value
            priorities = [pick_category(altitude_categories, altitude) for altitude in altitudes]
            segment_table["Altitude Value"] = altitudes
            segment_table["Altitude Priority"] = priorities


def calculate_list_priority(pl: PriorityList, list_name: str, invert: bool = False, false_value: float = 0.2) -> None:
    for _, row in pl.target_list.target_list.iterrows():
        target = row["Target Name"]
        list_value = row[f"List {list_name}"]
        if list_value != list_value: # set NaN values to False
            list_value = False
        priority = 1 if (list_value ^ invert) else false_value
        for table in pl.tables[target]:
            table[f"{list_name} Priority"] = priority


def calculate_prior_observation_priority(pl: PriorityList, prior_observation_categories: CategoryTable) -> None:
    pass


def calculate_phase_priority(pl: PriorityList, phase_defs: list[ph.PhaseEventDef], phase_categories: dict[str, float]) -> None:
    ephem_table = pl.target_list.other_lists["Ephem"]
    min_priority = min(phase_categories, key=operator.itemgetter(1))
    for target, table_list in pl.tables.items():
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
    for table_list in pl.tables.values():
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
    for target, target_tables in pl.tables.items():
        for target_table, aggregate_table in zip(target_tables, aggregate_tables):
            if np.max(target_table["Overall Priority"]) >= skip_column_threshold:
                aggregate_table[target] = target_table["Overall Priority"]
    final_tables = []
    for aggregate_table in aggregate_tables:
        aggregate_table.index = aggregate_table.index.strftime("%Y-%m-%d %H:%M")
        target_list = pl.target_list.target_list
        this_table_target_list = target_list[target_list["Target Name"].isin(aggregate_table.columns)]
        ra_order = list(this_table_target_list.sort_values("ra")["Target Name"])
        # ra_order = {this_table_target_list.sort_values("ra")[["Target Name", "ra"]]}
        final_tables.append(aggregate_table[ra_order])

    pl.target_priorities = final_tables


def aggregate_table_to_html(agg_table: pd.DataFrame, file: str = None) -> str:
    table_name = "targetList"
    caption = f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}"
    html = itables.to_html_datatable(
        df=agg_table,
        caption=caption,
        table_id=table_name,
        connected=True,
        paging=False,
        maxBytes=0,
        maxColumns=0,
        autoWidth=False,
        layout={"topStart": None, "topEnd": None},
        classes="display nowrap compact cell-border",
    )
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

    if file:
        with open(file, "w") as f:
            f.write(html)
    return html
