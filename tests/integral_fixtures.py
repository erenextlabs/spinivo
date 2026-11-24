"""Deterministic FIDs and analytic references for the integral evidence set."""

import json
from dataclasses import replace
from math import atan, pi
from pathlib import Path

import numpy as np

from spinivo.processing import ProcessingSettings
from spinivo.readers.model import RawFID

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SETTINGS = ProcessingSettings(zero_fill_factor=2, lb_hz=1, autophase_method="lineshape", smooth_baseline=True)
OBSERVE_MHZ, CARRIER_PPM = 400.0, 5.0
SW_HZ, POINTS = 4096, 32768
DEFAULT_BOUNDS_HZ = (-100, 100)
NOISE_BOUNDS_HZ = (300, 650)
DEFAULT_SNR = 3000
SEEDS = (1729, 1730, 1731, 1732, 1733, 1734, 1735, 1736)
TIME = np.arange(POINTS) / SW_HZ

# Each line is (center in Hz, relative whole-line area, intrinsic FWHM in Hz).
TARGETS = {
    "singlet": ((0, 1.0, 1),),
    "triplet": ((-7, 0.25, 1), (0, 0.5, 1), (7, 0.25, 1)),
    "satellites": ((-65, 0.0055, 1), (0, 0.989, 1), (65, 0.0055, 1)),
}
ANCHORS = ((-900, 1.0, 1), (1000, 1.0, 1))
FOREIGN = {"line": ((28, 0.05, 1),), "tail": ((110, 0.10, 29),)}


def load_cases():
    return json.loads((FIXTURES / "integral_cases.json").read_text(encoding="utf-8"))


def ppm_bounds(bounds_hz):
    return tuple(CARRIER_PPM + hz / OBSERVE_MHZ for hz in bounds_hz)


def processing_regions(case):
    return (ppm_bounds(case["bounds"]),)


def component_area(lines, bounds_hz):
    """Signed absorption area of ideal components after apodization, in intensity·Hz."""
    low, high = bounds_hz
    total = 0.0
    for center, area, width in lines:
        half_width = (width + SETTINGS.lb_hz) / 2
        total += area * (atan((high - center) / half_width) - atan((low - center) / half_width)) / pi
    return total


def target_lines(case):
    return TARGETS[case.get("target", "singlet")]


def declared_evidence(case):
    """Return the noise window and any extra-line assignment available to a labeller."""
    evidence = {"noise_window_hz": NOISE_BOUNDS_HZ}
    if "foreign" in case:
        evidence["extra_line_assignment"] = FOREIGN[case["foreign"]]
    if "withheld" in case:
        del evidence[case["withheld"]]
    return evidence


def _fid(lines):
    return sum(2 * area * np.exp(pi * (-width + 2j * center) * TIME) for center, area, width in lines)


def make_raw(case, seed=SEEDS[0], *, ideal=False):
    """Build a RawFID, or its noise-free and uncorrupted target-plus-anchors control."""
    target = target_lines(case)
    target_fid = _fid(target)
    fid = target_fid + case.get("anchor_scale", 1) * _fid(ANCHORS)
    if not ideal:
        if "phase_deg" in case:
            phase = np.exp(1j * np.deg2rad(case["phase_deg"]))
            fid += (phase - 1) * target_fid
        if "foreign" in case:
            fid += _fid(FOREIGN[case["foreign"]])
        if "baseline_area" in case:
            # Short Gaussian transient gives a broad, smooth absorption hump in the region.
            fid += 2 * case["baseline_area"] * np.exp(-2 * (pi * 20 * TIME) ** 2 + 2j * pi * 30 * TIME)
        weights = np.exp(-pi * SETTINGS.lb_hz * TIME)
        weights[0] = 0.5
        height = max(2 * area / (pi * (width + SETTINGS.lb_hz)) for _, area, width in target)
        noise_scale = height / case.get("snr", DEFAULT_SNR) * SW_HZ / np.linalg.norm(weights)
        rng = np.random.default_rng(seed)
        fid += noise_scale * (rng.standard_normal(POINTS) + 1j * rng.standard_normal(POINTS))
    return RawFID(
        vendor="synthetic", source=Path(case["id"]), nucleus="1H", observe_frequency_mhz=OBSERVE_MHZ,
        spectral_width_hz=SW_HZ, carrier_ppm=CARRIER_PPM, complex_points=POINTS,
        acquisition_time_s=POINTS / SW_HZ, relaxation_delay_s=None, pulse_width_us=None,
        flip_angle_deg=None, pulse_90_us=None, scans=None, dummy_scans=None, receiver_gain=None,
        temperature_k=None, pulse_program=None, digital_filter=None, intensity_note="Analytic fixture intensity",
        fid=fid, source_fields={}, missing=("relaxation_delay_s", "pulse_width_us", "flip_angle_deg", "pulse_90_us",
                                           "scans", "dummy_scans", "receiver_gain", "temperature_k", "pulse_program"),
    )


def make_baseline_raw(curvature=0.0):
    """Encode a Gaussian on a known background in a 4096-point, zero-fill-one FID."""
    raw = make_raw({"id": "baseline"}, ideal=True)
    hz = SW_HZ / 2 - np.arange(int(SW_HZ))
    trace = np.exp(-0.5 * (hz / 2) ** 2) + 0.25 + curvature * (hz / 1024) ** 2
    fid = np.fft.fft(np.fft.ifftshift(trace))
    fid[0] *= 2  # Invert first-point half weighting.
    return replace(raw, fid=fid, complex_points=fid.size, acquisition_time_s=fid.size / SW_HZ)
