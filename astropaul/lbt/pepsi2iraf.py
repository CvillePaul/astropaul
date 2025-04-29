from astropy.io import fits
import numpy as np
from scipy.interpolate import interp1d

def pepsi2iraf(in_file:str, out_file:str) -> None:

    # Open PEPSI FITS binary table
    hdul = fits.open(in_file)
    table = hdul[1].data

    wavelength = table["Arg"]
    flux = table["Fun"]

    # Create a linear wavelength grid
    start = wavelength[0]
    stop = wavelength[-1]
    step = np.median(np.diff(wavelength))  # approximate linear spacing
    wavelength_linear = np.arange(start, stop, step)

    # Interpolate flux onto new grid
    interp_flux = interp1d(wavelength, flux, kind="linear", bounds_error=False, fill_value="extrapolate")
    flux_linear = interp_flux(wavelength_linear)

    # Create a new 1D FITS image
    hdu = fits.PrimaryHDU(flux_linear)

    # Add WCS keywords to header
    hdu.header["CRVAL1"] = wavelength_linear[0]  # Starting wavelength
    hdu.header["CDELT1"] = step  # Step size
    hdu.header["CRPIX1"] = 1  # Reference pixel
    hdu.header["CTYPE1"] = "LINEAR"

    hdu.writeto(out_file, overwrite=True)
