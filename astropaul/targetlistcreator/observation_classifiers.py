import pandas as pd

import astropaul.phase as ph

from .targetlistcreator import TargetList, verify_step_requirements

observations_table = "PEPSI Observations"
results_table = "PEPSI Observation Results"
goals_table = "PEPSI Observation Goals"

def pepsi_observation_classifier(
    tl: TargetList = None,
    phase_event_defs: list[ph.PhaseEventDef] = None,
    min_evaluation: float = 0,
    **kwargs,
) -> TargetList:
    """Classify all PEPSI observations for each target.
    Identifier is Cross Disperser combined with a |-separated list of states for each system
    Merit is simply the Evaluation value for that observation"""
    verify_step_requirements(tl, {observations_table, "Ephemerides"})
    if not phase_event_defs:
        raise ValueError("Phase Event Defs not provided")
    answer = tl.copy()
    labels = pd.DataFrame(columns=["Target Name", "Label", "Merit"])
    all_ephems = tl.other_lists["Ephemerides"]
    quads = answer.target_list[answer.target_list["Target Type"] == "QuadEB"]["Target Name"].to_list()
    for target_name, rows in answer.other_lists[observations_table].groupby("Target Name"):
        if not target_name in quads:
            continue # only handle quads for now
        ephem_rows = all_ephems[(all_ephems["Target Name"] == target_name) & (all_ephems["Member"] == "a")].sort_values("System")
        target_ephems = [ph.Ephemeris.from_dataframe_row(ephem_row) for _, ephem_row in ephem_rows.iterrows()]
        for _, row in rows.iterrows():
            observation_time = row["Mid JD"]
            cross_disperser = f'CD{row["Cross Disperser"][0]}'
            phases = "|".join(
                [
                    ph.PhaseEventList.calc_phase_events(ephem, phase_event_defs, observation_time, observation_time).events[0].type
                    for ephem in target_ephems
                ]
            )
            label = f"{phases}_{cross_disperser}"
            merit = 1 if row["Evaluation"] >= min_evaluation else 0
            if merit != merit:
                continue
            labels.loc[len(labels)] = [target_name, label, merit]
    answer.other_lists[results_table] = labels.sort_values(["Target Name", "Label"])
    return answer


def add_pepsi_goals(tl: TargetList, **kwargs) -> TargetList:
    """Create a table that for each target lists the desired labels and total merit factor"""
    answer = tl.copy()
    # make list of all possible labels
    cross_dispersers = ["CD6"]
    phase_combos = ["R|T", "T|R", "B|T", "T|B", "T|T"]
    labels = [f"{phases}_{cd}" for phases in phase_combos for cd in cross_dispersers]
    # for each target, assign a desired merit for each label
    quads = tl.target_list[tl.target_list["Target Type"] == "QuadEB"]["Target Name"].to_list()
    goals = {target_name: {label: 4 for label in labels} for target_name in quads}
    # reduce the desired merit for each observation already made
    for target_name, rows in answer.other_lists[results_table].groupby("Target Name"):
        if not target_name in quads:
            continue
        target_goals = goals[target_name]
        for _, label, merit in rows[["Label", "Merit"]].itertuples():
            if not label in target_goals:
                continue
            current_goal = target_goals.get(label, 0)
            if current_goal <= merit:
                del target_goals[label]
            else:
                target_goals[label] = current_goal - merit
    rows = [[target_name, label, merit] for target_name, target_labels in goals.items() for label, merit in target_labels.items()]
    answer.other_lists[goals_table] = pd.DataFrame(rows, columns=["Target Name", "Label", "Merit"])
    return answer


# Original design notes
# * observation classifier function
#   * take a target list, an input table name and output table name
#   * input table is usually a table of observations
#   * output table indicates classifications, and has columns:
#     * target name
#     * identifier, which indicates a type of observation, such as:
#       * system configuration, like R|T indicating radial & transverse motion of a quad
#       * resolution, such as if it was observed with a 4m or 8m class scope
#       * spectral range, ie cross disperser used
#       * Some combination of these factors, eg CD2_R|T
#     * merit, a float that could indicate observation quality, integration time on target, or some other figure of merit
