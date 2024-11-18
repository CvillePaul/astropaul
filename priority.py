from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, get_body, GCRS
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd

import targetlistcreator as tlc
import phase as ph


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
        self.tables = {}
        for target in self.targets:
            segments = []
            for segment in self.segments:
                table = pd.DataFrame(index=[beg.to_datetime() for beg, _ in segment])
                segments.append(table)
            self.tables[target] = segments



CategoryTable = list[tuple[tuple[float, float], Any]]


def pick_category(categories: CategoryTable, value: float) -> Any:
    for (beg, end), result in categories:
        if beg <= value < end:
            return result


def calculate_moon_priority(
    pl: PriorityList, illumination_categories: CategoryTable = None, dist_categories: CategoryTable = None
) -> None:
    if not illumination_categories:
        illumination_categories = [
            ((0.0, 0.4), "Dark"),
            ((0.4, 0.7), "Gray"),
            ((0.7, 1.0), "Bright"),
        ]
    if not dist_categories:
        dist_categories = {
            "Dark": [
                ((0, 1), 1),
            ],
            "Gray": [
                ((0, 5), 0.1),
                ((5, 15), 0.5),
                ((15, 180), 1),
            ],
            "Bright": [
                ((0, 15), 0.1),
                ((15, 30), 0.5),
                ((30, 180), 1),
            ],
        }
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


def calculate_altitude_priority(pl: PriorityList, altitude_categories: CategoryTable = None) -> None:
    if not altitude_categories:
        altitude_categories = [
            ((-90, 35), 0),
            ((35, 45), 0.1),
            ((45, 60), 0.5),
            ((60, 90), 1),
        ]
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
        priority = (row[f"List {list_name}"] ^ invert) * (1 - false_value) + false_value
        for table in pl.tables[target]:
            table[f"{list_name} Priority"] = priority

def calculate_prior_observation_priority(pl: PriorityList, prior_observation_categories: CategoryTable) -> None:
    pass

def calculate_phase_priority(pl: PriorityList, phase_defs: list[ph.PhaseEventDef], ) -> None:
    pass

def calculate_overall_priority(pl: PriorityList, weightings: dict[str, float]) -> None:
    norm_factor = np.sum([weight for weight in weightings.values()])
    for table_list in pl.tables.values():
        for table in table_list:
            table["Overall Priority"] = 0
            for col, weight in weightings.items():
                table["Overall Priority"] += table[f"{col} Priority"] * weight
            table["Overall Priority"] /= norm_factor

