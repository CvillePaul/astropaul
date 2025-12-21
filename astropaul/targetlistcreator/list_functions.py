import collections
import inspect
from functools import partial

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import pandas as pd

from astropaul.database import string_to_db_style
from astropaul.observing import ObservingSession
import astropaul.phase as ph

from .targetlistcreator import TargetList, convert_columns_to_human_style, verify_step_requirements

import __main__

def add_targets(tl: TargetList = None, **kwargs) -> TargetList:
    """Add rows for each target in the database"""
    targets = pd.read_sql(
        "select * from targets;",
        kwargs["connection"],
    )
    convert_columns_to_human_style(targets)
    important_columns = ["Target Name", "RA", "Dec"]
    column_groups = {"Target": (important_columns, [column for column in targets.columns if not column in important_columns])}
    answer = TargetList(name=tl.name, target_list=targets, column_groups=column_groups)
    return answer


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
    list_columns = [f"List {name}" for name in target_lists]
    answer.column_groups["List"] = ([], list_columns)
    # also make a new table indicating all lists to which each target belongs
    lists_df = answer.target_list[["Target Name"] + list_columns]
    lists_df = lists_df.melt(id_vars="Target Name", var_name="List", value_name="value")
    lists_df = lists_df[lists_df["value"]].drop(columns="value")
    lists_df = lists_df.sort_values(["Target Name", "List"])
    answer.other_lists["List Memberships"] = lists_df
    return answer


def add_phase_events(
    tl: TargetList,
    observing_session: ObservingSession,
    phase_event_defs: list[ph.PhaseEventDef],
    event_types: list[str] = None,
    **kwargs,
) -> TargetList:
    """Calculate what phase events are in effect for each observing segment"""
    verify_step_requirements(tl, {"Ephemerides"})
    ephem_table = tl.other_lists["Ephemerides"]
    answer = tl.copy()
    ephems = {}
    phase_events = pd.DataFrame(
        columns=[
            "Target Name",
            "System",
            "Member",
            "Orbit Num",
            "Beg JD",
            "Phase",
            "Event",
            "End JD",
            "Beg UTC",
            "End UTC",
        ]
    )
    for _, ephem_row in ephem_table.sort_values(["System", "Member"]).iterrows():
        target_name = ephem_row["Target Name"]
        ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
        ephems[target_name] = ephems.get(target_name, []) + [ephem]
    for beg, end in observing_session.observing_segments:
        for target_name, ephem_list in ephems.items():
            for ephem in ephem_list:
                event_list = ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, beg, end).events
                event_list[0].time = beg  # treat first event as if it had started at segment beginning
                # set end time equal to begin time of next event, or if last event to segment end
                num_events = len(event_list)
                end_jds = [0.0] * num_events
                end_utcs = [""] * num_events
                for i in range(num_events):
                    if i == num_events - 1:
                        this_end = end  # make sure we don't go past end of observing segment
                    else:
                        this_end = event_list[i + 1].time
                    end_jds[i] = this_end
                    end_utcs[i] = Time(this_end, format="jd").iso[:19]

                # populate table with all values except end times
                for event, this_end, end_utc in zip(event_list, end_jds, end_utcs):
                    if event_types and not event.type in event_types:
                        continue  # skip unwanted phase events
                    phase_events.loc[len(phase_events)] = (
                        [target_name] + event.to_list() + [this_end, event.time.iso[:19] if event.time else "", end_utc]
                    )
    answer.other_lists["Phase Events"] = phase_events.sort_values(["Target Name", "System", "Member", "Beg JD"])
    return answer


def ancillary_data_from_tess(tl: TargetList, **kwargs) -> TargetList:
    """Add fundamental things like magnitude and effective temperature, if available, from TESS catalog"""
    verify_step_requirements(tl)
    answer = tl.copy()
    tess_data = pd.read_sql("select * from tess_ticv8;", kwargs["connection"])
    if len(tess_data) > 0:
        convert_columns_to_human_style(tess_data)
        primary_columns = ["Vmag", "Teff"]
        answer.target_list = answer.target_list.merge(
            tess_data[["Target Name"] + primary_columns], on="Target Name", how="left"
        )
        answer.other_lists["TESS"] = tess_data
        answer.column_groups["TESS Data"] = (primary_columns, [])
    return answer


def add_columns_from_sql(tl: TargetList, table_name: str, primary_cols: list[str], **kwargs) -> TargetList:
    """Join columns from another table to the target list.  Other table must have a Target Name column."""
    verify_step_requirements(tl)
    answer = tl.copy()
    other_table = pd.read_sql(f"select * from {string_to_db_style(table_name)};", kwargs["connection"])
    if not other_table.empty:
        convert_columns_to_human_style(other_table)
        if not "Target Name" in other_table.columns:
            raise ValueError(f"Table {table_name} must have a column named Target Name")
        secondary_cols = [col for col in other_table.columns if col not in primary_cols + ["Target Name"]]
        answer.target_list = answer.target_list.merge(other_table, on="Target Name", how="left")
        answer.column_groups[table_name] = (primary_cols, secondary_cols)
    return answer


def hide_cols(tl: TargetList, prefix: str, **kwargs) -> TargetList:
    answer = tl.copy()
    cols_to_remove = [col for col in tl.target_list.columns if col.startswith(prefix)]
    if len(cols_to_remove) > 0:
        answer.target_list.drop(labels=cols_to_remove, axis=1, inplace=True)
    del answer.column_groups[prefix]
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
        code = code[prefix + 9 : -2]  # try to isolate the code of the lambda function passed in
    # replace local variables used in the lambda function with their values
    for variable in criteria.__code__.co_names:
        if value := __main__.__dict__.get(variable, None):
            code = code.replace(variable, repr(value))
    if inverse:
        code = "Inverse of: " + code
        answer.target_list = answer.target_list[~criteria(answer.target_list)]
    else:
        answer.target_list = answer.target_list[criteria(answer.target_list)]
    answer.list_criteria.add(code)
    return answer


def filter_other_list(tl: TargetList, list_name: str, criteria, inverse: bool = False, **kwargs) -> TargetList:
    verify_step_requirements(tl, {list_name})
    answer = tl.copy()
    other_list = Table.from_pandas(tl.other_lists[list_name])
    if len(other_list) == 0:
        return answer
    code = inspect.getsource(criteria).strip()
    if (prefix := code.find("criteria=")) > 0:
        code = f"In table {list_name}: {code[prefix + 9 :-2]}"  # try to isolate the code of the lambda function passed in
    # replace local variables used in the lambda function with their values
    for variable in criteria.__code__.co_names:
        if variable == "u":
            continue  # hacky way to prevent the u alias for astropy.units from looking ugly in output
        if value := __main__.__dict__.get(variable, None):
            code = code.replace(variable, repr(value))
    # convert to astropy table, with each column containing a Quantity having its unit attribute set
    # all so the dear user can use lambdas like df["col"] > 1.2 * u.arcsec, which isn't possible on dataframes
    # NOTE: this requires all rows to have EXACTLY the same unit
    for col in other_list.columns:
        sample_val = other_list[col][0]
        if isinstance(sample_val, u.Quantity):
            unit = sample_val.unit
            other_list[col] = [val.value for val in other_list[col]]
            other_list[col].unit = unit
    mask = criteria(other_list)
    if inverse:
        code = "Inverse of: " + code
        mask = ~mask
    surviving_targets = set(other_list[mask]["Target Name"])
    answer.target_list = answer.target_list[tl.target_list["Target Name"].isin(surviving_targets)]
    answer.list_criteria.add(code)
    return answer


def tabulate_ephemerides(tl: TargetList, **kwargs) -> TargetList:
    """Redo counting of ephemeris rows: 1.0 for rows with all values specified, 0.1 for rows lacking a duration"""
    verify_step_requirements(tl, {"Ephemerides"})
    answer = tl.copy()
    ephem = tl.other_lists["Ephemerides"].copy()
    ephem["score"] = ephem["Duration"].apply(lambda x: 1.0 if x == x else 0.1)
    ephem_scores = ephem.groupby("Target Name", as_index=False)["score"].sum()
    ephem_scores.columns = ["Target Name", "Num Ephemerides"]
    answer.target_list.drop("Num Ephemerides", axis=1, inplace=True)
    answer.target_list = answer.target_list.merge(ephem_scores, on="Target Name", how="left")
    answer.target_list["Num Ephemerides"] = answer.target_list["Num Ephemerides"].fillna(0.0)
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


def add_member_phases(
    tl: TargetList,
    table_name: str,
    time_column: str,
    synthetic_phase_percent: float = float("nan"),
    **kwargs,
) -> TargetList:
    """For the specified table, add a column indicating the phases of each member of the target.
    For these calculations, the specified column contains astropy Time objects to use.
    A phase of 0.0 is mid eclipse of the a member.
    If `synthetic_phase_percent` is provided it will be used on systems lacking an eclipse duration."""
    verify_step_requirements(tl, {table_name})
    table = tl.other_lists[table_name]
    for col in ["Target Name", time_column]:
        if not col in table.columns:
            raise ValueError(f"Table {table_name} does not have a column {time_column}")
    answer = tl.copy()
    all_targets = set(table["Target Name"].unique())
    ephem = tl.other_lists["Ephemerides"]
    ephem = ephem[ephem["Target Name"].isin(all_targets)]
    ephem = ephem.sort_values(["Target Name", "System", "Member"])
    target_ephems = collections.defaultdict(list)
    for _, row in ephem.iterrows():
        target_name = row["Target Name"]
        target_ephems[target_name].append(ph.Ephemeris.from_dataframe_row(row))

    phases = []
    for _, row in table.iterrows():
        target_name = row["Target Name"]
        observation_time = row[time_column]
        phases.append(
            " ".join(
                [
                    f"{target_ephem.system}{target_ephem.component}: {ph.calc_phase(target_ephem, observation_time):.2f}"
                    for target_ephem in target_ephems[target_name]
                ]
            )
        )
    table["Member Phases"] = phases
    answer.other_lists[table_name] = table
    return answer


def add_system_configuration(
    tl: TargetList,
    table_name: str,
    time_column: str,
    synthetic_phase_percent: float = float("nan"),
    eclipse_table: str = None,
    **kwargs,
) -> TargetList:
    """For the specified table, add a pair of columns indicating the state of each system of the target.
    For these calculations, the specified column contains astropy Time objects to use.
    For a system in eclipse, show what member and the percentage of the time between ingress and egress.
    For a system not in eclipse, show the percent of phase, where 0.0 is mid eclipse of the a member.
    If `synthetic_phase_percent` is provided it will be used on systems lacking an eclipse duration.
    If `eclipse_table` is specified, add a table of all observations occurring during eclipse."""
    verify_step_requirements(tl, {table_name, "Ephemerides"})
    table = tl.other_lists[table_name]
    for col in ["Target Name", time_column]:
        if not col in table.columns:
            raise ValueError(f"Table {table_name} does not have a column {time_column}")
    answer = tl.copy()
    all_targets = set(table["Target Name"].unique())
    # determine max number of systems we need to handle, max duration
    ephem = tl.other_lists["Ephemerides"]
    ephem = ephem[ephem["Target Name"].isin(all_targets)]
    systems = sorted(ephem["System"].unique())
    # add the desired columns to the table
    phase_event_defs = [
        ph.PhaseEventDef(
            "Not in Eclipse", partial(ph.calc_time_of_gress, ingress=False, synthetic_phase_percent=synthetic_phase_percent)
        ),
        ph.PhaseEventDef(
            "Eclipse", partial(ph.calc_time_of_gress, ingress=True, synthetic_phase_percent=synthetic_phase_percent)
        ),
    ]
    # create a separate table of observations made during eclipse to add to answer if requested
    eclipse_observations = pd.DataFrame(columns=["Target Name", "Mid JD", "Mid UTC", "System", "Member", "Type"])
    # categorize all observation rows
    system_eclipses, system_phases = {}, {}
    for system in systems:
        system_eclipses[system] = []
        system_phases[system] = []
    for _, row in table.iterrows():
        target_name = row["Target Name"]
        observation_time = row[time_column]
        target_ephem = ephem[ephem["Target Name"] == target_name]
        for system in systems:
            phases = ""
            system_ephem = target_ephem[target_ephem["System"] == system]
            if system_ephem.empty:  # this target doesn't have an ephem entry for this system
                system_eclipses[system].append("")
                continue
            eclipse_found = False
            for _, ephem_row in system_ephem.sort_values("Member").iterrows():
                member_ephem = ph.Ephemeris.from_dataframe_row(ephem_row)
                phases += f"{system}{ephem_row["Member"]}: {ph.calc_phase(member_ephem, observation_time):.2f} "
                if member_ephem.duration != member_ephem.duration:
                    continue  # skip calculation for members with no duration specified
                event_list = ph.PhaseEventList.calc_phase_events(
                    member_ephem, phase_event_defs, observation_time, observation_time + member_ephem.period
                )
                if len(event_list.events) < 2:
                    raise ValueError(f"Unexpected event list {event_list}")
                if event_list.events[0].type == "Eclipse":
                    duration_percent = 1 - ((event_list.events[1].time - observation_time) / member_ephem.duration).to(
                        u.dimensionless_unscaled
                    )
                    system_eclipses[system].append(f"{system}{ephem_row["Member"]}: {duration_percent:.2f}")

                    eclipse_found = True
                    eclipse_observations.loc[len(eclipse_observations)] = [
                        target_name,
                        observation_time.jd,
                        observation_time.iso[:19],
                        ephem_row["System"],
                        ephem_row["Member"],
                        "Real duration" if member_ephem.duration == member_ephem.duration else "Synthetic duration",
                    ]
                    break
            if not eclipse_found:
                system_eclipses[system].append("")
            system_phases[system].append(phases)
        for system in systems:
            table[f"System {system} Eclipse"] = system_eclipses[system]
            table[f"System {system} Phases"] = system_phases[system]
    answer.other_lists[table_name] = table
    if eclipse_table:
        answer.other_lists[eclipse_table] = eclipse_observations
        # add a count column to target list indicating how many rows per target in eclipse_table
        count_column = f"Num {eclipse_table}"
        counts = eclipse_observations["Target Name"].value_counts()
        answer.target_list[count_column] = answer.target_list["Target Name"].map(counts).fillna(0).astype(int)
        primary_columns, secondary_columns = answer.column_groups.get("Count", ([], []))
        primary_columns.append(count_column)
        answer.column_groups["Count"] = (primary_columns, secondary_columns)
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
    for idx, pepsi in pepsi_table.iterrows():
        ephem_rows = ephem_table[ephem_table["Target Name"] == pepsi["Target Name"]]
        if ephem_rows.empty:
            continue
        exposure_days = pepsi["Exposure Time"]
        beg = Time(pepsi["Mid JD"], format="jd") - exposure_days / 2
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
        pepsi_table.loc[idx, "RV Status"] = "|".join(system_states)
        # pepsi_phases.loc[len(pepsi_phases)] = [
        #     pepsi["Target Name"],
        #     beg,
        #     Time(beg, format="jd").iso[:18],
        #     "|".join(system_states),
        # ]
    # if not pepsi_phases.empty:
    #     answer.other_lists["PEPSI RV Status"] = pepsi_phases
    answer.other_lists["PEPSI Observations"] = pepsi_table
    return answer


def add_pepsi_evaluations(tl: TargetList, **kwargs) -> TargetList:
    table_name = "PEPSI Observations"
    verify_step_requirements(tl, {table_name})
    answer = tl.copy()
    evaluations = pd.read_sql("select * from pepsi_evaluations;", kwargs["connection"])
    if len(evaluations) > 0:
        convert_columns_to_human_style(evaluations)
        pepsi_observations = answer.other_lists[table_name]
        pepsi_observations = pepsi_observations.merge(
            evaluations[["Spectrum File", "Evaluation"]], on="Spectrum File", how="left"
        )
        answer.other_lists[table_name] = pepsi_observations
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


def add_tess_catalog_associations(tl: TargetList, **kwargs) -> TargetList:
    # use several of the catalog associations available in the TESS results table
    answer = tl.copy()
    conn = kwargs["connection"]
    all_targets = list(tl.target_list["Target Name"].unique())
    catalog_members = pd.read_sql("select target_name, catalog_name, catalog_id from catalog_members;", conn)
    tess_catalog_ids = pd.read_sql(
        """
        select cm.target_name, tt.gaia 'GAIA DR2', hip, twomass
        from catalog_members cm
        join tess_ticv8 tt on tt.tic = cm.catalog_id
        where cm.catalog_name = 'TESS TICv8';
   
    """,
        conn,
    )
    tess_members = tess_catalog_ids.melt(id_vars="target_name", var_name="catalog_name", value_name="catalog_id")
    tess_members = tess_members[tess_members["catalog_id"] == tess_members["catalog_id"]]
    catalog_members = pd.concat([catalog_members, tess_members], ignore_index=True)
    convert_columns_to_human_style(catalog_members)
    catalog_members = catalog_members[catalog_members["Target Name"].isin(all_targets)]
    answer.other_lists["Catalog Membership"] = catalog_members
    return answer


def add_database_table(tl: TargetList, table_name: str, add_count: bool = True, **kwargs) -> TargetList:
    verify_step_requirements(tl)
    answer = tl.copy()
    conn = kwargs["connection"]
    sql_name = string_to_db_style(table_name)
    table_contents = pd.read_sql(f"select * from {sql_name};", conn)
    # apply any transformations, such as units
    table_metadata = pd.read_sql(
        f"select column_name, value from table_metadata where table_name = '{sql_name}' and value_type='unit';", conn
    )
    for column_name, unit in table_metadata[["column_name", "value"]].values:
        match unit:
            case "JD":
                table_contents[column_name] = [Time(jd, format="jd") for jd in table_contents[column_name]]
            case _:  # by default, assume the unit is an astropy unit
                table_contents[column_name] = [u.Quantity(value, unit=unit) for value in table_contents[column_name]]
    convert_columns_to_human_style(table_contents)
    # for tables that indicate target name, restrict table to targets in this TargetList
    if "Target Name" in table_contents.columns:
        all_targets = list(tl.target_list["Target Name"].unique())
        table_contents = table_contents[table_contents["Target Name"].isin(all_targets)]
        if add_count:  # optionally add a column to the TargetList indicating num rows (in this table) per target
            count_column = f"Num {table_name}"
            counts = table_contents["Target Name"].value_counts()
            answer.target_list[count_column] = answer.target_list["Target Name"].map(counts).fillna(0).astype(int)
            primary_columns, secondary_columns = answer.column_groups.get("Count", ([], []))
            primary_columns.append(count_column)
            answer.column_groups["Count"] = (primary_columns, secondary_columns)
    answer.other_lists[table_name] = table_contents
    return answer


def add_resource(tl: TargetList, resource: str, table_name: str, **kwargs) -> TargetList:
    verify_step_requirements(tl, {table_name})
    answer = tl.copy()
    conn = kwargs["connection"]
    resources = pd.read_sql(f"select * from resource_{string_to_db_style(resource)};", conn)
    convert_columns_to_human_style(resources)
    answer.other_lists[table_name] = answer.other_lists[table_name].merge(resources, on="ID", how="left")
    return answer
