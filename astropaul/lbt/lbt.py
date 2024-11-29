import pandas as pd

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
        max(
            (times := pepsi_exptime(vmag, snr, teff=teff, fiber_setup=fiber, binocular=binocular))[cd_red - 1],
            times[cd_blue - 1],
            60,  # never recommend exposures under 1 minute
        )
        for vmag, teff in answer.target_list[["Vmag", "Teff"]].values
    ]
    answer.target_list[f"{column_prefix}priority"] = priority
    answer.target_list[f"{column_prefix}notes"] = notes
    primary_columns = {f"{column_prefix}exp_time", f"{column_prefix}priority", f"{column_prefix}notes"}
    secondary_columns = set([col for col in answer.target_list.columns if col.startswith(column_prefix)]) - primary_columns
    answer.column_groups[column_prefix] = (list(primary_columns), list(secondary_columns))
    return answer


def make_lbt_readme_table(target_list: pd.DataFrame) -> pd.DataFrame:
    targets = target_list.copy().sort_values("ra")
    # TODO: put code to pivot rows to a given lst
    readme = pd.DataFrame()
    readme["Target Name"] = targets["Target Name"]
    readme["RA"] = targets["RA"]
    readme["Dec"] = targets["Dec"]
    readme["Vmag"] = targets["Vmag"]
    readme["Teff"] = targets["Teff"]
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


def write_lbt_readme_file(file_base: str, targets: pd.DataFrame) -> str:
    readme = make_lbt_readme_table(targets)
    # save target list as csv
    readme.to_csv(
        file_base + ".csv",
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
