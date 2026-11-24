"""Measure declared signal and noise regions on a processed spectrum."""

from dataclasses import dataclass
from math import isfinite

import numpy as np

from .acquisition import SignalMeasurement
from .processing import BaselineRegion, Spectrum, integrate


@dataclass(frozen=True)
class RegionMeasurement:
    """Signed area, peak and noise measurements with units and availability reasons."""

    name: str
    signal_bounds_ppm: tuple[float, float]
    noise_bounds_ppm: tuple[float, float] | None
    area: float | None
    peak_ppm: float | None
    height: float | None
    noise_sd: float | None
    fwhm_hz: float | None
    source: str
    transformed_points: int
    area_reason: str | None = None
    peak_ppm_reason: str | None = None
    height_reason: str | None = None
    noise_sd_reason: str | None = None
    fwhm_reason: str | None = None
    area_unit: str = "intensity*Hz"
    peak_unit: str = "ppm"
    height_unit: str = "intensity"
    noise_sd_unit: str = "intensity"
    fwhm_unit: str = "Hz"
    baseline_region: BaselineRegion | None = None

    def to_signal_measurement(self) -> SignalMeasurement:
        """Return the values and provenance used by the acquisition checks."""
        return SignalMeasurement(
            height=self.height, noise_sd=self.noise_sd, fwhm_hz=self.fwhm_hz, source=self.source,
            noise_bounds_ppm=self.noise_bounds_ppm, transformed_points=self.transformed_points,
            fwhm_reason=self.fwhm_reason,
        )


def _noise_sd(ppm, real, signal_bounds, noise_bounds):
    if noise_bounds is None:
        return None, "Missing noise window"
    low, high = noise_bounds
    if not isfinite(low) or not isfinite(high):
        return None, "Noise bounds must be finite"
    signal_low, signal_high = signal_bounds
    if signal_low <= signal_high and signal_low <= high and low <= signal_high:
        return None, "Noise window overlaps the signal bounds"
    start = ppm.searchsorted(low)
    stop = ppm.searchsorted(high, side="right")
    if stop - start < 2:
        return None, "Noise window contains fewer than two spectrum points"
    if low < ppm[0] or high > ppm[-1]:
        return None, "Noise window reaches outside the computed spectral grid"
    return float(np.std(real[start:stop], ddof=1)), None


def _half_height_crossing(hz, real, half_height):
    for index, value in enumerate(real[1:], start=1):
        if value <= half_height:
            before = real[index - 1]
            fraction = (before - half_height) / (before - value)
            return hz[index - 1] + fraction * (hz[index] - hz[index - 1])
    return None


def _fwhm(spectrum, start, stop, peak, bounds, peak_reason):
    if spectrum.settings.lb_hz > 0:
        return None, "FWHM requires lb_hz = 0"
    if peak is None:
        return None, peak_reason
    ppm, hz, real = spectrum.ppm[::-1], spectrum.hz[::-1], spectrum.data.real[::-1]
    half_height = real[peak] / 2
    # Include the adjacent samples beyond the bounds to bracket crossings between grid points.
    left = slice(max(0, start - 1), peak + 1)
    right = slice(peak, stop + 1)
    low_crossing = _half_height_crossing(hz[left][::-1], real[left][::-1], half_height)
    high_crossing = _half_height_crossing(hz[right], real[right], half_height)
    low_hz, high_hz = np.interp(bounds, ppm, hz)
    if (low_crossing is None or high_crossing is None
            or low_crossing < low_hz or high_crossing > high_hz):
        return None, "Missing a half-height crossing inside the signal bounds"
    return float(high_crossing - low_crossing), None


def measure_region(spectrum: Spectrum, name: str, signal_bounds_ppm: tuple[float, float],
                   noise_bounds_ppm: tuple[float, float] | None) -> RegionMeasurement:
    """Measure a declared isolated line and an optional signal-free window on the real spectrum."""
    low, high = sorted(signal_bounds_ppm)
    signal_bounds = (low, high)
    noise_bounds = None
    if noise_bounds_ppm is not None:
        noise_low, noise_high = sorted(noise_bounds_ppm)
        noise_bounds = (noise_low, noise_high)
    ppm, real = spectrum.ppm[::-1], spectrum.data.real[::-1]

    area, area_reason = None, None
    try:
        area = integrate(spectrum, signal_bounds).area
    except ValueError as error:
        area_reason = str(error)

    start = ppm.searchsorted(low)
    stop = ppm.searchsorted(high, side="right")
    outside_grid = low < ppm[0] or high > ppm[-1]
    peak, peak_ppm, height, peak_reason = None, None, None, None
    if outside_grid:
        peak_reason = "Signal window reaches outside the computed spectral grid"
    elif not isfinite(low) or not isfinite(high) or low >= high:
        peak_reason = "Signal bounds must span a finite nonzero interval"
    elif start == stop:
        peak_reason = "Signal window contains no spectrum points"
    else:
        candidate = start + np.argmax(real[start:stop])
        if real[candidate] > 0:
            peak, peak_ppm, height = candidate, float(ppm[candidate]), float(real[candidate])
        else:
            peak_reason = "Signal window maximum must be positive"

    noise_sd, noise_reason = _noise_sd(ppm, real, signal_bounds, noise_bounds)
    if outside_grid:
        fwhm, fwhm_reason = None, peak_reason
    else:
        fwhm, fwhm_reason = _fwhm(spectrum, start, stop, peak, signal_bounds, peak_reason)
    baseline = spectrum.smooth_baseline_fit
    diagnostic = None
    if baseline is not None:
        diagnostic = next((region for region in baseline.regions if region.bounds_ppm == signal_bounds), None)
    return RegionMeasurement(
        name=name, signal_bounds_ppm=signal_bounds, noise_bounds_ppm=noise_bounds,
        area=area, peak_ppm=peak_ppm, height=height, noise_sd=noise_sd, fwhm_hz=fwhm,
        source=f"{spectrum.source}: processed spectrum ({name})", transformed_points=spectrum.transformed_points,
        area_reason=area_reason, peak_ppm_reason=peak_reason, height_reason=peak_reason,
        noise_sd_reason=noise_reason, fwhm_reason=fwhm_reason, baseline_region=diagnostic,
    )
