import astroplan as ap
from astropy.coordinates import SkyCoord, get_body
import astropy.units as u
import numpy as np
import pandas as pd

from astropaul.observing import ObservingSession

from .targetlistcreator import TargetList

def add_observability(
    tl: TargetList,
    observing_session: ObservingSession,
    observability_threshold: tuple[u.Quantity, u.Quantity] = (30 * u.deg, 80 * u.deg),
    constraints: list[ap.Constraint] = None,
    column_prefix: str = "Observable ",
    calc_moon_distance: bool = False,
    time_resolution=15 * u.min,
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
    location = observing_session.observer.name or observing_session.observer.location
    time_range = f"{observing_session.time_range[0].iso[:10]} to {observing_session.time_range[1].iso[:10]}"
    answer.list_criteria.add(f"Observability calculated at {location} in {time_resolution} intervals from {time_range}")
    for constraint in constraints:
        answer.list_criteria.add(f"{constraint.__class__.__name__}: {vars(constraint)}")
    overall_any_night = np.array([False] * len(answer.target_list))
    overall_every_night = np.array([True] * len(answer.target_list))
    overall_max_alts = np.array([-90.0] * len(answer.target_list))
    overall_min_hours_observable = np.array([24] * len(answer.target_list))
    overall_min_moon_dist = np.array([180.0] * len(answer.target_list))
    nightly_columns = []
    new_cols = {}
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
        new_cols[column_name] = segment_observability
        overall_any_night |= segment_observability
        overall_every_night &= segment_observability
        nightly_columns.append(column_name)
        hours_observable_column = f"{column_name} Hours Observable"
        overall_min_hours_observable = [
            np.min([current_min_hours, this_min_hours])
            for current_min_hours, this_min_hours in zip(overall_min_hours_observable, segment_hours_observable)
        ]
        new_cols[hours_observable_column] = segment_hours_observable
        nightly_columns.append(hours_observable_column)
        max_alt_column = f"{column_name} Max Alt"
        new_cols[max_alt_column] = segment_max_alts
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
            new_cols[moon_dist_column] = dist
            nightly_columns.append(moon_dist_column)
            overall_min_moon_dist = [
                np.min([current_min_dist, this_dist]) for current_min_dist, this_dist in zip(overall_min_moon_dist, dist)
            ]
    any_night_column = f"{column_prefix}Any Night"
    new_cols[any_night_column] = overall_any_night
    every_night_column = f"{column_prefix}Every Night"
    new_cols[every_night_column] = overall_every_night
    main_columns = [any_night_column, every_night_column]
    max_alt_column = f"{column_prefix}Max Alt"
    new_cols[max_alt_column] = overall_max_alts
    main_columns.append(max_alt_column)
    min_hours_observable_column = f"{column_prefix}Hours Min"
    new_cols[min_hours_observable_column] = overall_min_hours_observable
    main_columns.append(min_hours_observable_column)
    if calc_moon_distance:
        moon_dist_column = f"{column_prefix}Min Moon Dist"
        new_cols[moon_dist_column] = overall_min_moon_dist
        main_columns.append(moon_dist_column)
    answer.target_list = answer.target_list.assign(**new_cols)
    answer.column_groups[column_prefix.strip()] = (main_columns, nightly_columns)
    moon_phases = pd.DataFrame()
    moon_phases["Time"] = [beg.iso[:10] for beg, _ in observing_session.observing_segments]
    moon_phases["Phase"] = [ap.moon_illumination(t) for t, _ in observing_session.observing_segments]
    answer.add_other("Lunar Phases", moon_phases)
    return answer
