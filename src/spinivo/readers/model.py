"""Common result and error types for raw FID readers."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class RawFID:
    """Raw time-domain samples and acquisition metadata in common units.

    intensity_scale is the multiplier the reader applied to the stored values; None means unknown.
    storage_scale identifies the vendor's storage scaling convention, including the Bruker NC exponent.
    """

    vendor: str
    source: Path
    nucleus: str
    observe_frequency_mhz: float | None
    spectral_width_hz: float | None
    carrier_ppm: float | None
    complex_points: int
    acquisition_time_s: float | None
    relaxation_delay_s: float | None
    pulse_width_us: float | None
    flip_angle_deg: float | None
    pulse_90_us: float | None
    scans: int | None
    dummy_scans: int | None
    receiver_gain: float | None
    temperature_k: float | None
    pulse_program: str | None
    digital_filter: dict | None
    intensity_note: str
    fid: NDArray[np.complex128]
    source_fields: dict[str, str]
    missing: tuple[str, ...]
    intensity_scale: float | None = None
    storage_scale: str | None = None


class UnsupportedInput(ValueError):
    """An input that cannot be read as a supported raw FID."""

    def __init__(self, reason, source, detail):
        self.reason = reason
        self.source = source
        self.detail = detail
        super().__init__(f"{source}: {detail}")


def _validate_fields(source, params, *, integers=(), reals=(), strings=()):
    for keys, types in ((integers, (int,)), (reals, (int, float)), (strings, (str,))):
        for key in keys:
            value = params.get(key)
            if value is None:
                continue
            try:
                valid = type(value) in types and (float not in types or isfinite(value))
            except OverflowError:
                valid = False
            if not valid:
                raise UnsupportedInput("truncated", source, f"invalid value for {key}")


def _require_proton(source, key, nucleus, proton_names):
    name = "" if nucleus is None else nucleus.strip().casefold()
    if not name:
        raise UnsupportedInput("missing_nucleus", source, f"missing {key}")
    if name not in proton_names:
        raise UnsupportedInput("not_proton", source, f"{key} is {nucleus!r}")


def _raw_fid(vendor, source, fid, intensity_note, metadata, *, intensity_scale, storage_scale):
    values = {field: None if isinstance(value, str) and not value.strip() else value
              for field, (value, _) in metadata.items()}
    sources = {field: key for field, (_, key) in metadata.items() if values[field] is not None}
    width = values["spectral_width_hz"]
    values["acquisition_time_s"] = None
    if width is not None:
        if not isfinite(width) or width <= 0:
            raise UnsupportedInput("truncated", source, f'{sources["spectral_width_hz"]} must be finite and positive')
        values["acquisition_time_s"] = values["complex_points"] / width
        if not isfinite(values["acquisition_time_s"]):
            raise UnsupportedInput("truncated", source, "acquisition time is out of range")
        sources["acquisition_time_s"] = f'{sources["complex_points"]}, {sources["spectral_width_hz"]}'
    return RawFID(
        vendor=vendor, source=source, intensity_note=intensity_note, intensity_scale=intensity_scale,
        storage_scale=storage_scale,
        fid=np.asarray(fid, dtype=np.complex128), source_fields=sources,
        missing=tuple(field for field, value in values.items() if value is None), **values,
    )
