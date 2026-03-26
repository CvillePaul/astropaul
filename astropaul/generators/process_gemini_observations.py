from io import StringIO
from pathlib import Path

import pandas as pd  

# json files were downloaded from the archive search page at https://archive.gemini.edu/searchform/defaults/ by program ID

# Our program IDs:
# 2024B: GN-2024B-Q-305 and GS-2024B-Q-308
# 2025A: GN-2025A-Q-307 and GS-2025A-Q-305
# 2025B: proposal not approved (DARP rules)
# 2026A: only GS-2026A-Q-307 (no observations from Gemini North)

def process_gemini_observations(observations_json: str) -> pd.DataFrame:
    sequences = pd.read_json(StringIO(observations_json))
    columns_to_keep = [
        "telescope",
        "instrument",
        "ut_datetime",
        "object",
        "ra",
        "dec",
        "airmass",
        "exposure_time",
    ]
    sequences = sequences[columns_to_keep].sort_values(["instrument", "ut_datetime"])
    sequences["ut_date"] = [date[:11] for date in sequences["ut_datetime"]]
    observations = sequences[["telescope", "instrument", "object", "ut_date"]].drop_duplicates()
    # observations = observations[observations["object"].str.startswith("TIC")]
    # sorted(observations["ut_date"].unique())
    observations.columns = ["Telescope", "Instrument", "Target Name", "Date UTC"]
    return observations

def process_all_gemini_observations(input_dir: str, output_dir: str) -> None:
    for file in Path(input_dir).glob("*.json"):
        with open(file, "r") as f:
            contents = f.read()
        observations = process_gemini_observations(contents)
        outfile = Path(output_dir) / f"{file.name}.csv"
        observations.to_csv(outfile, index=False)
        print(f"Wrote {outfile}")
