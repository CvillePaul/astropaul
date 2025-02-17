import collections
from typing import Any

import astroplan as ap
from astropy.coordinates import SkyCoord, Angle, get_body
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd

from astropaul.observing import ObservingSession
import astropaul.phase as ph


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

    def add_columns(self, columns: dict[str, Any] = {}):
        pass

    def add_other(self, name: str, other: pd.DataFrame):
        self.other_lists[name] = other

    def summarize(self) -> str:
        answer = f"{self.name}\n"
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
    count_column = f"Num {column_prefix.replace(" ", "")}"
    num_speckle = pd.read_sql(
        f"""
        select ss.target_id, count(ss.id) as '{count_column}'
        from tom_specklesession ss
        group by ss.target_id
        ;
        """,
        kwargs["connection"],
        index_col="target_id",
    )
    column_groups = {"Count": ([count_column], [])}
    # next, add a separate table of all speckle observations
    speckle = pd.read_sql(
        """
        select t.name 'Target Name', ss.speckle_session Session, ss.num_sequences Sequences, ss.start_jd Start, ss.mid_jd Mid, ss.end_jd End
        from tom_specklesession ss
        join tom_target as t on t.id = ss.target_id
        ;
        """,
        kwargs["connection"],
        index_col="Target Name",
    )
    existing_targets = tl.target_list["Target Name"].to_list()
    speckle = speckle[speckle.index.isin(existing_targets)]
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in speckle.columns}
    speckle.rename(columns=new_column_names, inplace=True)
    speckle[f"{column_prefix}Mid UTC"] = Time(speckle[f"{column_prefix}Mid"], format="jd").iso
    answer = TargetList.merge(tl, num_speckle, column_groups, {"Speckle": speckle})
    answer.target_list[count_column] = answer.target_list[count_column].fillna(int(0))
    return answer


def add_pepsi(tl: TargetList, column_prefix="PEPSI ", **kwargs) -> TargetList:
    # first, add a column for total observations to main table
    count_column = f"Num {column_prefix.replace(" ", "")}"
    spectra_count = pd.read_sql(
        f"""
        select target_id, count(id) as '{count_column}'
        from tom_spectrumrawdata
        group by target_id
        ;""",
        kwargs["connection"],
        index_col="target_id",
    )
    column_groups = {"Count": ([count_column], [])}
    # next, add a separate table of all spectra
    spectra = pd.read_sql(
        """
        select t.name 'Target Name', srd.*
        from tom_spectrumrawdata as srd
        join tom_target as t on t.id = srd.target_id
        """,
        kwargs["connection"],
        index_col="Target Name",
    )
    existing_targets = tl.target_list["Target Name"].to_list()
    spectra = spectra[spectra.index.isin(existing_targets)]
    new_column_names = {column: f"{column_prefix}{column.capitalize()}" for column in spectra.columns}
    spectra.rename(columns=new_column_names, inplace=True)
    answer = TargetList.merge(tl, spectra_count, column_groups, {"PEPSI": spectra})
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


def add_ephemerides(tl: TargetList, **kwargs) -> TargetList:
    ephem = pd.read_sql(
        """
        select t.name "Target Name", bp.system, bp.member, bp.period, bp.t0, bp.duration, bp.depth
        from tom_binaryparameters bp
        join tom_scienceresult sr on sr.id = bp.scienceresult_ptr_id
        join tom_target t on t.id = sr.target_id
        ;""",
        kwargs["connection"],
        index_col="Target Name",
    )
    existing_targets = tl.target_list["Target Name"]
    ephem = ephem[ephem.index.isin(existing_targets)]
    ephem["period"] = ephem["period"] / 3600 / 24  # convert from seconds to days for convenience
    ephem["duration"] = ephem["duration"] / 3600 / 24  # convert from seconds to days for convenience
    new_column_names = {column: column.capitalize() for column in ephem.columns}
    ephem.rename(columns=new_column_names, inplace=True)
    answer = tl.copy()
    answer.add_other("Ephem", ephem)
    return answer


def add_phase_events(
    tl: TargetList, observing_session: ObservingSession, phase_event_defs: list[ph.PhaseEventDef], **kwargs
) -> TargetList:
    """Calculate what phase events are in effect for each observing segment"""
    if (ephem_table := tl.other_lists.get("Ephem", pd.DataFrame())).empty:
        raise ValueError("Ephem table not present")
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
                event_list[0].jd = beg.jd # treat first event as if it had started at segment beginning
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
                        [target_name]
                        + event.to_list()
                        + [end_jd, Time(event.jd, format="jd").iso[:19], end_utc]
                    )
    answer.other_lists["Phase Events"] = phase_events
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
    coords = SkyCoord(ra=tl.target_list["ra"].values * u.deg, dec=tl.target_list["dec"].values * u.deg)

    # calculate observability for each time segment of the observing session
    answer = tl.copy()
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


def add_speckle_phase(tl: TargetList, phase_event_defs: list[ph.PhaseEventDef], **kwargs) -> TargetList:
    """Calculate what PhaseEventDef was most in effect during a speckle observation"""
    required_tables = {"Speckle", "Ephem"}
    existing_tables = set(tl.other_lists.keys())
    if not existing_tables.issuperset(required_tables):
        raise ValueError(f"One or more required table missing: {', '.join(required_tables - existing_tables)}")
    answer = tl.copy()
    ephem_table = tl.other_lists["Ephem"]
    speckle_phases = pd.DataFrame(columns=["Target Name", "Speckle Session", "JD Mid", "UTC Mid", "System", "Member", "State"])
    for target_name, speckle in tl.other_lists["Speckle"].iterrows():
        ephem_rows = ephem_table[ephem_table.index == target_name]
        if ephem_rows.empty:
            continue
        beg, end = speckle["Speckle Start"], speckle["Speckle End"]
        for _, ephem_row in ephem_rows.sort_values(["System", "Member"]).iterrows():
            ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
            state = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg, end).calc_longest_span(beg, end)
            speckle_phases.loc[len(speckle_phases)] = [
                target_name,
                int(speckle["Speckle Session"]),
                speckle["Speckle Mid"],
                Time(speckle["Speckle Mid"], format="jd").iso[:18],
                ephem_row["System"],
                ephem_row["Member"],
                state,
            ]
    if not speckle_phases.empty:
        answer.other_lists["Speckle Phases"] = speckle_phases
    return answer
