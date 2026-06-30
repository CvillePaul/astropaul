from io import StringIO
from pathlib import Path

import pandas as pd  

# json files were downloaded from the archive search page at https://archive.gemini.edu/searchform/defaults/ by program ID

# Our program IDs:
# 2024B: GN-2024B-Q-305 and GS-2024B-Q-308
# 2025A: GN-2025A-Q-307 and GS-2025A-Q-305
# 2025B: proposal not approved (DARP rules)
# 2026A: only GS-2026A-Q-307 (no observations from Gemini North)
# 2026B: GN-2026B-Q-308 and GS-2026B-Q-308

def process_gemini_observation_file(observation_file: Path, output_dir: Path, verbose: bool = False) -> pd.DataFrame:
    with open(observation_file, "r") as f:
        observations_json = f.read()
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
    csv_file = output_dir / f'{observation_file.name.replace(".json", ".csv")}'
    observations.to_csv(csv_file, index=False)
    if verbose:
        print(f"Wrote {len(observations)} lines to {csv_file}")
