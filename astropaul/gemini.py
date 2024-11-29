import astropy.units as u
import pandas as pd

def add_gemini_speckle_params(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    # For more info, see https://www.gemini.edu/instrumentation/alopeke-zorro/proposal-preparation
    iq70_exposures, iq85_exposures = [], []
    for vmag in df["Vmag"].values:
        if vmag < 9:
            iq70_exposures.append(3)
            iq85_exposures.append(6)
        elif 9 <= vmag < 12:
            iq70_exposures.append(5)
            iq85_exposures.append(10)
        elif 12 <= vmag < 14:
            iq70_exposures.append(7)
            iq85_exposures.append(14)
        elif 14 <= vmag < 16:
            iq70_exposures.append(11)
            iq85_exposures.append(22)
        elif 16 <= vmag:
            iq70_exposures.append(15)
            iq85_exposures.append(30)
        elif vmag != vmag: # test for nan
            iq70_exposures.append(15)
            iq85_exposures.append(30)
        else:
            raise ValueError(f"Cannot handle Vmag of {vmag}")
    df["gemini_iq70_speckle_sequences"] = iq70_exposures
    df["gemini_iq85_speckle_sequences"] = iq85_exposures
    speckle_sequence_duration = 64 * u.s
    df["gemini_iq70_speckle_exposure"] = df["gemini_iq70_speckle_sequences"] * speckle_sequence_duration
    df["gemini_iq85_speckle_exposure"] = df["gemini_iq85_speckle_sequences"] * speckle_sequence_duration
    pointing_time = 300 * u.s
    calibaration_overhead = 1.167 # 10 minutes of calibration target per hour of science target
    photometric_overhead = 1.016 # 10 minutes of calibration per 10 hour observing session
    total_overhead = calibaration_overhead * photometric_overhead
    df["gemini_iq70_speckle_program_time"] = (df["gemini_iq70_speckle_exposure"] + pointing_time)* total_overhead
    df["gemini_iq85_speckle_program_time"] = (df["gemini_iq85_speckle_exposure"] + pointing_time)* total_overhead
    return df
