from datetime import datetime
import operator
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, get_body
from astropy.time import Time
import astropy.units as u
import itables
import numpy as np
import pandas as pd

from astropaul.html import dataframe_to_datatable
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
        self.categorized_priorities: list[pd.DataFrame] = None
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
        self.categorized_priorities = []
        for numerical_priority in self.numerical_priorities:
            category_table = numerical_priority.apply(
                lambda col: pd.cut(col, bins=bins, labels=labels, ordered=False, include_lowest=True)
            )
            self.categorized_priorities.append(category_table)

    def to_html(self, file_base: str = "Priorities", dir: str = "."):
        segment_starts = [segment[0][0].iso[:10] for segment in self.segments]
        table_options = {"search": False, "sort": False}
        # make a page for each target's priority scores
        for target, target_table_list in self.target_tables.items():
            for start_utc, target_table in zip(segment_starts, target_table_list):
                tt = target_table.copy()
                tt.index = [f"{time:%H:%M}" for time in tt.index]
                target_html = f'<h1 style="text-align: center">{start_utc} Target Scores for {target}</h1>\n'
                target_html += dataframe_to_datatable(tt, "targetPriority", table_options=table_options)
                with open(f"{dir}/Target {file_base} {target} {start_utc}.html", "w") as f:
                    f.write(target_html)
        # make a page for numerical priorities for each observing segment
        if self.numerical_priorities:
            for start_utc, priority_table in zip(segment_starts, self.numerical_priorities):
                pt = priority_table.copy()
                pt.index = [f"{time:%H:%M}" for time in pt.index]
                threshold = self.category_bins[-2]
                for col in pt.columns:
                    pt[col] = [
                        val if val < threshold else f'<p style="background-color: #EBF4FA">{val:.3f}</p>' for val in pt[col]
                    ]
                pt.columns = [f'<a href="Target {file_base} {target} {start_utc}.html">{target}</a>' for target in pt.columns]
                numerical_html = f'<h1 style="text-align: center">{start_utc} Target Priorities</h1>\n'
                numerical_html += dataframe_to_datatable(pt, "numericalPriority", table_options=table_options)
                with open(f"{dir}/Numerical {file_base} {start_utc}.html", "w") as f:
                    f.write(numerical_html)
        # make a page for categorized priorities for each observing segment
        if self.categorized_priorities:
            for start_utc, categories_table in zip(segment_starts, self.categorized_priorities):
                ct = categories_table.copy()  # make a copy we can alter for formatting purposes
                ct.index = [f"{time:%H:%M}" for time in ct.index]
                list_targets = self.target_list.target_list
                # not all targets on list have nonzero priority, so we need a list of targets for this particular segment
                segment_targets = list_targets[list_targets["Target Name"].isin(ct.columns)]
                # add some helpful informational rows to the top of the chart
                ct.loc["RA"] = [
                    segment_targets[segment_targets["Target Name"] == col]["RA"].values[0][:8] for col in ct.columns
                ]
                ct.loc["Dec"] = [
                    segment_targets[segment_targets["Target Name"] == col]["Dec"].values[0][:9] for col in ct.columns
                ]
                ct.loc["Teff"] = [
                    f"{segment_targets[segment_targets["Target Name"] == col]["Teff"].values[0]:.0f}" for col in ct.columns
                ]
                new_ordering = [*ct.index[-3:], *ct.index[:-3]]
                ct = ct.reindex(new_ordering)
                # elide the target names to first 4 digits
                ct.columns = [col[:8] + "..." for col in ct.columns]
                # generate html version of the chart
                categorized_html = f'<h1 style="text-align: center">{start_utc} Target Priorities</h1>\n'
                categorized_html += dataframe_to_datatable(ct, "categorizedPriority", table_options=table_options)
                with open(f"{dir}/Categorical {file_base} {start_utc}.html", "w") as f:
                    f.write(categorized_html)
        # make summary page
        summary_html = f"<h1>{self.target_list.name}</h1>\n"
        summary_html += '<table cellpadding="10" border="1">\n'
        summary_html += (
            f"<tr><td>Observatory</td><td>{self.session.observer.name}, {self.session.observer.timezone}</td></tr>\n"
        )
        summary_html += (
            f'<tr><td>Target List</td><td><a href="{self.target_list.name}.html">{self.target_list.name}</a></td></tr>\n'
        )
        summary_html += (
            f'<tr><td>Session Start</td><td style="font-family: monospace">{self.session.time_range[0].iso}</td></tr><br>\n'
        )
        summary_html += (
            f'<tr><td>Session Finish</td><td style="font-family: monospace">{self.session.time_range[1].iso}</td></tr><br>\n'
        )
        summary_html += f"<tr><td>Table Interval</td><td>{self.interval}</td></tr>\n"
        summary_html += "</table><br>\n"
        summary_html += "<h2>Target List Summary</h2>\n"
        summary_html += '<pre style="font-family: monospace">\n'
        summary_html += self.target_list.summarize()
        summary_html += "</pre>\n"
        summary_html += "<br><br><h2>Observing Segments</h2>\n"
        summary_html += '<table cellpadding="10" style="text-align: center;">'
        summary_html += "    <tr><th>Start</th><th>Finish</th><th>Numerical Priorities</th><th>Categorical Priorities</th></tr>"
        for segment, segment_start in zip(self.segments, segment_starts):
            summary_html += "<tr>"
            summary_html += f"<td>{segment[0][0].iso[:19]}</td><td>{segment[-1][1].iso[:19]}</td>"
            summary_html += f'<td><a href="Numerical {file_base} {segment_start}.html">Link</a></td>'
            summary_html += f'<td><a href="Categorical {file_base} {segment_start}.html">Link</a></td>'
            summary_html += "</tr>\n"
        summary_html += "</table>\n"
        with open(f"{dir}/Summary {self.target_list.name}.html", "w") as f:
            f.write(summary_html)


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
        first_ra_of_night = lst - 60  # lst is the RA at 90 deg alt, we want targets at 30 deg alt
        # now "pivot" the table so the first column is the RA
        i = next(i for i, ra in enumerate(this_target_list["ra"].tolist()) if ra > first_ra_of_night)  # 1st viewable target
        target_names = target_names[i:] + target_names[:i]  # "rotate" the list so first ra is one > lst
        pl.numerical_priorities.append(aggregate_table[target_names])
