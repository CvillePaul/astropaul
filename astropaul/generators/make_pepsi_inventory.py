import argparse
from glob import glob
import os
import sys

from astropy.io import fits
from astropy.table import Table

spectrum_extensions = [".nor", ".bwl"]

header_items = {
    "DATE-OBS": ("Date UTC", str),
    "TIME-OBS": ("Exposure Start UTC", str),
    "OBJECT": ("Target Name", str),
    "RA2000": ("RA", str),
    "DE2000": ("Dec", str),
    "TELESCOP": ("Observatory", str),
    "INSTRUME": ("Instrument", str),
    "ARM": ("Arm", str),
    "CROSDIS": ("Cross Disperser", str),
    "FIBER": ("Fiber", str),
    "EXPTIME": ("Exposure Time", float),
    "JD-OBS": ("Mid JD", float),
}


def make_pepsi_inventory(spectrum_dirs: list[str], out_dir: str = ".", verbose: bool = False):
    col_names, col_types = zip(*header_items.values())
    col_names = list(col_names)
    col_names.append("Spectrum File")
    col_types = list(col_types)
    col_types.append("str")
    for spectrum_dir_pattern in spectrum_dirs:
        fixed_dir_pattern = spectrum_dir_pattern.replace('"', '') # powershell does crappy things w/ trailing backslashes
        for spectrum_dir in glob(fixed_dir_pattern):
            pepsi_observations = Table(names=col_names, dtype=col_types)
            spectrum_files = sum(
                [glob(os.path.join(spectrum_dir, f"pepsi*{extension}")) for extension in spectrum_extensions], []
            )
            for spectrum_file in spectrum_files:
                try:
                    hdr = fits.open(spectrum_file)[0].header
                    values = [field_type(hdr[field_name]) for field_name, (_, field_type) in header_items.items()]
                    values.append(os.path.basename(spectrum_file))
                    pepsi_observations.add_row(values)
                    if verbose:
                        print(f"Processed file: {spectrum_file}")
                except Exception as e:
                    print(f"Exception while processing file: {spectrum_file}")
                    print(e)
                    return
            if len(pepsi_observations) > 0:
                observations_file = f"PEPSI Observations {os.path.basename(spectrum_dir)}.csv"
                pepsi_observations.write(os.path.join(out_dir, observations_file), overwrite=True)
                if verbose:
                    print(f"Wrote {len(pepsi_observations)} observations to {observations_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description="Use FITS headers to make a CSV file of PEPSI observations for each specified directory.",
    )
    parser.add_argument("-o", "--out_dir", default=".", help="Output directory (default = %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output stats about file processing")
    parser.add_argument("spectrum_dir", nargs="*", help="List of directories to process, wildcards allowed")
    args = parser.parse_args()
    make_pepsi_inventory(args.spectrum_dir, args.out_dir, verbose=args.verbose)
