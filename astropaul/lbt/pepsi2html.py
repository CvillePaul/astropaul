import argparse
from glob import glob
from datetime import datetime, timedelta
import os
from astropy.io import fits
from bokeh.io import output_file
from bokeh.layouts import column
from bokeh.plotting import figure, output_file, save
from bokeh.models.tools import BoxZoomTool, PanTool, ResetTool
import numpy as np

spectrum_extensions = [".nor", ".bwl"]

def pepsi2html(spectrum_dirs: list[str], out_dir: str = ".", verbose: bool = False) -> None:
    for spectrum_dir_pattern in spectrum_dirs:
        fixed_dir_pattern = spectrum_dir_pattern.replace('"', '') # powershell does crappy things w/ trailing backslashes
        for spectrum_dir in glob(fixed_dir_pattern):
            dir_name = os.path.basename(os.path.normpath(spectrum_dir))
            outfile = os.path.join(out_dir, f"PEPSI Spectra {dir_name}.html")
            output_file(outfile)
            plots = []
            spectrum_files = sum(
                [glob(os.path.join(spectrum_dir, f"pepsi*{extension}")) for extension in spectrum_extensions], []
            ) # the sum with [] is a trick to flatten a list of lists
            for spectrum_file in spectrum_files:
                filename = os.path.basename(spectrum_file)
                hdul = fits.open(spectrum_file)
                header = hdul[0].header
                target = header["OBJECT"]
                obstime = f"{header["DATE-OBS"]} {header["TIME-OBS"]}"
                arm = header["ARM"]
                cd = header["CROSDIS"]
                exposure = header["EXPTIME"]
                if not type(exposure) is float:
                    dt = datetime.strptime(exposure, "%H:%M:%S")
                    exposure = timedelta(hours=dt.hour, minutes=dt.minute, seconds=dt.second).total_seconds()
                if (snr := header.get("SNR")) is None:
                    snr = 0

                data = hdul[1].data
                title = (
                    f"{target:<23s}: {arm:>4s} arm, CD{cd:<14s}, {exposure:>4.0f}sec, SNR {snr:5d}, {obstime[:19]}  File: {filename:<30s}"
                )
                p = figure(
                    title=title,
                    sizing_mode="stretch_width",
                    tools=[BoxZoomTool(), PanTool(), ResetTool()],
                )
                p.line(data["Arg"], np.clip(data["Fun"], a_min=None, a_max=1.3), line_width=0.2)
                p.xaxis.ticker.num_minor_ticks = 25
                plots.append(p)
                # if verbose:
                #     print(title)

            save(column(plots, sizing_mode="scale_width"))
            # to convert to PDF, open html file in browser and print to file
            if verbose:
                print(f"Output saved to {outfile}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="pepsi2html",
        description="Convert reduced PEPSI .bwl file(s) into a single HTML file showing all spectra",
    )

    parser.add_argument("-o", "--outdir", default=".", help="Directory for html output file (default = %(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output file info as each file is processed")
    parser.add_argument("filenames", nargs="+", help="List of PEPSI spectrum files to process")
    args = parser.parse_args()
    pepsi2html(args.filenames, args.outdir, args.verbose)
