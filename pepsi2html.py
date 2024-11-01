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

parser = argparse.ArgumentParser(
    prog="pepsi2html",
    description="Convert reduced PEPSI .bwl file(s) into a single HTML file showing all spectra",
)

parser.add_argument("-o", "--outfile", default="pepsi_spectra.html", help="Name of file to create (default = %(default)s)")
parser.add_argument("-v", "--verbose", action="store_true", help="Output file info as each file is processed")
parser.add_argument("filenames", nargs="+", help="List of PEPSI spectrum files to process")
args = parser.parse_args()

output_file(args.outfile)
plots = []
for filepattern in args.filenames:
    for path in glob(filepattern):
        filename = os.path.basename(path)
        hdul = fits.open(path)
        header = hdul[0].header
        target = header["OBJECT"]
        obstime = f"{header["DATE-OBS"]} {header["TIME-OBS"]}"
        telescope = header["TELESCOP"]
        arm = header["ARM"]
        fiber = header["FIBER"]
        cd = header["CROSDIS"]
        exposure = header["EXPTIME"]
        if not type(exposure) is float:
            dt = datetime.strptime(exposure, "%H:%M:%S")
            exposure = timedelta(hours=dt.hour, minutes=dt.minute, seconds=dt.second).total_seconds()
        if (snr := header.get("SNR")) is None:
            snr = 0

        data = hdul[1].data
        title = (
            f"{target:<20s}: {arm:>4s} arm, CD{cd:<14s}, {exposure:>4.0f}sec, SNR {snr:5d}, {obstime[:19]}  File: {filename:<30s}"
        )
        p = figure(
            title=title,
            sizing_mode="stretch_width",
            tools=[BoxZoomTool(), PanTool(), ResetTool()],
        )
        p.line(data["Arg"], np.clip(data["Fun"], a_min=None, a_max=1.3), line_width=0.2)
        p.xaxis.ticker.num_minor_ticks = 25
        # p.add_layout(Label(x=np.min(data["Arg"]), y=1.2, x_units="data", y_units="data", text=f"Quality={quality}"))

        # output_file(filename=f"{target}_{arm}.html", title=p.title)
        # f = f"{target}_{arm[:1]}_{obstime[:13]}h{obstime[14:16]}m{obstime[17:19]}s.html"
        # output_file(filename=f, title=p.title)
        plots.append(p)
        # save(p)
        if args.verbose:
            print(title)

save(column(plots, sizing_mode="scale_width"))
# to convert to PDF, open html file in browser and print to file
if args.verbose:
    print(f"Output saved to {args.outfile}")
