from dataclasses import dataclass

from astropy.convolution import Box1DKernel, convolve
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from specutils import Spectrum
from specutils.analysis import correlation
from specutils.fitting import find_lines_derivative, find_lines_threshold
from specutils.manipulation import extract_region
from specutils.spectra import SpectralRegion


@dataclass
class PepsiSpectrum:
    spectrum: Spectrum
    target_name: str
    obs_time: Time
    coord: SkyCoord
    location: EarthLocation

    def copy(self) -> "PepsiSpectrum":
        return PepsiSpectrum(
            spectrum=Spectrum(self.spectrum.flux, self.spectrum.spectral_axis),
            target_name=self.target_name,
            obs_time=self.obs_time,
            coord=self.coord,
            location=self.location,
        )

    def extract_region(self, region: SpectralRegion):
        if region:
            self.spectrum = extract_region(self.spectrum, region)

    def smooth(self, boxcar: int = 25):
        self.spectrum = Spectrum(
            convolve(self.spectrum.flux, Box1DKernel(boxcar), boundary="extend"), self.spectrum.spectral_axis
        )

    def find_lines(self, flux_threshold: float, boxcar: float = 25):
        smoothed_spectrum = self.copy()
        smoothed_spectrum.smooth(boxcar=boxcar)
        return find_lines_derivative(smoothed_spectrum.spectrum, flux_threshold=flux_threshold)

    def correlate(self, other_spectrum: "PepsiSpectrum", region: SpectralRegion = None):
        sub_spectrum = extract_region(self.spectrum, region) if region else self.spectrum
        sub_other = extract_region(other_spectrum.spectrum, region) if region else other_spectrum.spectrum
        return correlation.template_correlate(sub_spectrum, sub_other)

    @staticmethod
    def read(filename: str, radial_velocity: u.Quantity = 0 * u.km / u.s) -> "PepsiSpectrum":
        hdul = fits.open(filename)
        header = hdul[0].header
        target_name = header["OBJECT"]
        location = EarthLocation(lat=header["LATITUDE"], lon=header["LONGITUD"], height=header["ALTITUDE"])
        obs_time = Time(f"{header["DATE-OBS"]} {header["TIME-OBS"]}", format="iso")
        coord = SkyCoord(ra=header["RA2000"], dec=header["DE2000"], unit=(u.hourangle, u.deg), frame="icrs")
        barycentric_rv = coord.radial_velocity_correction(kind="barycentric", obstime=obs_time, location=location)
        data = hdul[1].data
        wavelength, flux = data["Arg"] * u.AA, data["Fun"] * u.dimensionless_unscaled
        spectrum = Spectrum(spectral_axis=wavelength, flux=flux - 1, radial_velocity=barycentric_rv)
        spectrum.shift_spectrum_to(radial_velocity=spectrum.radial_velocity + radial_velocity)
        return PepsiSpectrum(spectrum=spectrum, target_name=target_name, obs_time=obs_time, coord=coord, location=location)
