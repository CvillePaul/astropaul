from pathlib import Path

from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoLocator, AutoMinorLocator
import pandas as pd
from specutils import Spectrum

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
    def find_standard(target_coord: SkyCoord, teff: float, standards: pd.DataFrame, standards_coords: SkyCoord) -> int:
        nearby_standards = standards[target_coord.separation(standards_coords) < 60 * u.deg]
        print(len(nearby_standards))
        closest_idx = (nearby_standards["Teff"] - teff).abs().idxmin()
        if closest_idx == closest_idx:
            return nearby_standards.loc[closest_idx, "Target Name"]
        else:
            return ""

    answer = tl.copy()
    targets = answer.target_list
    rv_standards = targets[targets["Target Type"] == "RV Standard"]
    standards_coords = SkyCoord(ra=rv_standards["RA"], dec=rv_standards["Dec"], unit="deg")

    targets["RV Standard"] = [
        (
            find_standard(SkyCoord(ra=ra, dec=dec, unit="deg"), teff, rv_standards, standards_coords)
            if target_type in target_types
            else ""
        )
        for target_type, ra, dec, teff in targets[["Target Type", "RA", "Dec", "Teff"]].values
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


def read_pepsi_file(filename: str) -> Spectrum:
    hdul = fits.open(filename)
    header = hdul[0].header
    target_name = header["OBJECT"]
    location = EarthLocation(lat=header["LATITUDE"], lon=header["LONGITUD"], height=header["ALTITUDE"])
    obs_time = Time(f"{header["DATE-OBS"]} {header["TIME-OBS"]}", format="iso")
    coord = SkyCoord(ra=header["RA2000"], dec=header["DE2000"], unit=(u.hourangle, u.deg), frame="icrs")
    barycentric_rv = coord.radial_velocity_correction(kind="barycentric", obstime=obs_time, location=location)
    # barycentric_rv = header["SSBVEL"] * u.m / u.s
    data = hdul[1].data
    wavelength, flux = data["Arg"] * u.AA, data["Fun"] * u.dimensionless_unscaled
    spectrum = Spectrum(spectral_axis=wavelength, flux=flux - 1, meta={"header": header})
    spectrum.shift_spectrum_to(radial_velocity=barycentric_rv)
    spectrum.meta["Target Name"] = target_name
    spectrum.meta["Location"] = location
    spectrum.meta["Observation Time"] = obs_time
    spectrum.meta["Coord"] = coord
    spectrum.meta["Barycenter RV"] = barycentric_rv
    spectrum.meta["Filename"] = filename
    return spectrum

def plot_pepsi_spectrum(spectrum: Spectrum):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(spectrum.spectral_axis, spectrum.flux, linewidth=0.1)
    cross_disperser = spectrum.meta["header"]["CROSDIS"]
    xmin = int(cross_disperser[3:7])
    xmax = int(cross_disperser[10:14])
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=9, prune=None))
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=10))
    ax.set_ylabel("Flux (au)")
    ax.margins(y=0.15)
    ax.grid(axis="x", color="0.70", which="major")
    ax.grid(axis="x", color="0.95", which="minor")
    filename = Path(spectrum.meta["Filename"]).name
    ax.set_title(f"{spectrum.meta["Target Name"]}   {spectrum.meta["Observation Time"].iso[:19]}   {filename}")
    return fig, ax
