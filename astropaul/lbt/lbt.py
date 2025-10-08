from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd

import astropaul.observing as obs
import astropaul.targetlistcreator as tlc
from .exposure_time import pepsi_exptime


def add_pepsi_params(
    tl: tlc.TargetList,
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


def assign_rv_standards(tl: tlc.TargetList, target_types: set[str], drop_unused: bool = True, **kwargs) -> pd.DataFrame:
    def find_standard(value: float, standards: pd.DataFrame) -> int:
        closest_idx = (standards["Teff"] - value).abs().idxmin()
        if closest_idx == closest_idx:
            return standards.loc[closest_idx, "Target Name"]
        else:
            return ""

    answer = tl.copy()
    targets = answer.target_list
    rv_standards = targets[targets["Target Type"] == "RV Standard"]

    targets["RV Standard"] = [
        find_standard(teff, standards=rv_standards) if target_type in target_types else ""
        for target_type, teff in targets[["Target Type", "Teff"]].values
    ]
    if drop_unused:
        all_standards = set(rv_standards["Target Name"])
        standards_to_drop = all_standards - set(targets["RV Standard"])
        targets = targets[~targets["Target Name"].isin(standards_to_drop)]
    answer.target_list = targets
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
    coords = SkyCoord(ra=targets["RA"], dec=targets["Dec"], unit=u.deg)
    # build up the output table
    readme = pd.DataFrame()
    readme["Target Name"] = targets["Target Name"]
    readme["RA"] = coords.ra.to_string(unit=u.hour, sep=":", precision=2)
    readme["Dec"] = coords.dec.to_string(unit=u.hour, sep=":", precision=2, alwayssign=True)
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
