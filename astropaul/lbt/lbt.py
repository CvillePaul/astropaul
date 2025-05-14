import pandas as pd

import astropaul.observing as obs
from astropaul.targetlistcreator import TargetList
from .exposure_time import pepsi_exptime


def add_pepsi_params(
    tl: TargetList,
    fiber: str,
    cd_blue: int,
    cd_red: int,
    snr: int,
    binocular: bool = False,
    priority: str = "",
    notes: str = "",
    column_prefix="PEPSI ",
    **kwargs,
) -> pd.DataFrame:
    answer = tl.copy()
    answer.target_list[f"{column_prefix}fiber"] = fiber
    answer.target_list[f"{column_prefix}cd_blue"] = cd_blue
    answer.target_list[f"{column_prefix}cd_blue_num_exp"] = 1
    answer.target_list[f"{column_prefix}cd_red"] = cd_red
    answer.target_list[f"{column_prefix}cd_red_num_exp"] = 1
    answer.target_list[f"{column_prefix}snr"] = snr
    answer.target_list[f"{column_prefix}exp_time"] = [
        round(
            max(
                (times := pepsi_exptime(vmag, snr, teff=teff, fiber_setup=fiber, binocular=binocular))[cd_red - 1],
                times[cd_blue - 1],
                60,  # never recommend exposures under 1 minute
            )
        )
        for vmag, teff in answer.target_list[["Vmag", "Teff"]].values
    ]
    answer.target_list[f"{column_prefix}priority"] = priority
    answer.target_list[f"{column_prefix}notes"] = notes
    primary_columns = {f"{column_prefix}exp_time", f"{column_prefix}priority", f"{column_prefix}notes"}
    secondary_columns = set([col for col in answer.target_list.columns if col.startswith(column_prefix)]) - primary_columns
    answer.column_groups[column_prefix] = (list(primary_columns), list(secondary_columns))
    return answer


def assign_rv_standards(
    tl: TargetList, science_types: list[str], rv_standard_type: str = "RV Standard", **kwargs
) -> pd.DataFrame:
    answer = tl.copy()
    rv_standards = answer.target_list[answer.target_list["Target Type"] == rv_standard_type]
    targets = answer.target_list
    # science_targets = answer.target_list[answer.target_list["Target Type"].isin(science_types)]

    def find_standard(value: float, standards: pd.DataFrame) -> int:
        closest_idx = (standards["Teff"] - value).abs().idxmin()
        if closest_idx == closest_idx:
            return standards.loc[closest_idx, "Target Name"]
        else:
            return ""

    targets["RV Standard"] = [
        find_standard(teff, standards=rv_standards) if type in science_types else ""
        for teff, type in targets[["Teff", "Target Type"]].values
    ]
    # science_targets["Teff"].apply(find_standard, standards=rv_standards)
    # answer.target_list = answer.target_list.merge(science_targets[["Teff", "RV Standard"]], on="Teff", how="left")
    # answer.target_list["RV Standard"].fillna("", inplace=True)
    return answer


def make_lbt_readme_table(target_list: pd.DataFrame, beg_lst: float = 0) -> pd.DataFrame:
    targets = target_list.copy().sort_values("RA")
    # "pivot" the rows of the table so it starts with the first ra > beg_lst
    i = 0
    for _, row in targets.iterrows():
        if row["RA"] > beg_lst:
            break
        i += 1
    targets = pd.concat([targets.iloc[i:], targets.iloc[:i]], ignore_index=True)
    # build up the output table
    readme = pd.DataFrame()
    readme["Target Name"] = targets["Target Name"]
    readme["RA"] = targets["RA HMS"]
    readme["Dec"] = targets["Dec DMS"]
    readme["Vmag"] = targets["Vmag"]
    readme["Teff"] = [f"{val:.0f}" for val in targets["Teff"]]
    readme["Fiber"] = targets["PEPSI fiber"]
    readme["BLUE Cross Disperser"] = targets["PEPSI cd_blue"]
    readme["BLUE CD NExp"] = targets["PEPSI cd_blue_num_exp"]
    readme["BLUE CD Exp Time"] = [f"{val:.0f}" for val in targets["PEPSI exp_time"]]
    readme["RED Cross Disperser"] = targets["PEPSI cd_red"]
    readme["RED CD NExp"] = targets["PEPSI cd_red_num_exp"]
    readme["RED CD Exp Time"] = [f"{val:.0f}" for val in targets["PEPSI exp_time"]]
    readme["Desired BLUE SNR"] = targets["PEPSI snr"]
    readme["Desired RED SNR"] = targets["PEPSI snr"]
    readme["Priority"] = targets["PEPSI priority"]
    readme["Notes"] = targets["PEPSI notes"]
    return readme


def write_lbt_readme_file(file_base: str, targets: pd.DataFrame, session: obs.ObservingSession) -> str:
    readme = make_lbt_readme_table(targets, session.starting_lst - 60)
    # save target list as csv
    readme.to_csv(
        file_base + ".csv",
        index=False,
    )
    # make the readme file by prepending/appending the header/footer info
    readme_header = open(file_base + ".README.header", "r").readlines()
    readme_footer = open(file_base + ".README.footer", "r").readlines()
    output = ""
    for line in readme_header:
        output += line.rstrip() + "\n"
    output += readme.to_string(index=False)
    for line in readme_footer:
        output += line.rstrip() + "\n"
    with open(file_base + ".README", "w") as f:
        f.write(output)
    return output
