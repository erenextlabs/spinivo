"""Resolve and remove group delay from raw 1D proton FIDs."""

from dataclasses import dataclass
from math import ceil, isfinite, prod

import numpy as np
from nmrglue.fileio.bruker import bruker_dsp_table
from numpy.typing import NDArray

from .readers.model import RawFID


@dataclass(frozen=True)
class DelayDecision:
    """Group delay in acquired complex points and the metadata used."""

    delay_points: float | None
    method: str
    source_fields: tuple[tuple[str, object], ...]
    reason: str | None = None


@dataclass(frozen=True)
class DelayRemoval:
    """Corrected samples and the applied delay, tail trim and first-point treatment."""

    fid: NDArray[np.complex128]
    decision: DelayDecision
    method: str
    tail_trim_points: int
    first_point_treatment: str


def _source_fields(metadata, *keys):
    fields = []
    for key in keys:
        value = metadata.get(key)
        fields.append((key, tuple(value) if isinstance(value, list) else value))
    return tuple(fields)


def _positive_integer(value):
    return type(value) is int and value > 0


def _bruker_delay(metadata):
    grpdly = metadata.get("GRPDLY")
    try:
        valid = type(grpdly) in (float, int) and isfinite(grpdly) and grpdly > 0
    except OverflowError:
        valid = False
    if valid:
        return DelayDecision(float(grpdly), "GRPDLY", _source_fields(metadata, "GRPDLY"))

    source = _source_fields(metadata, "GRPDLY", "DSPFVS", "DECIM")
    dspfvs, decim = metadata.get("DSPFVS"), metadata.get("DECIM")
    if _positive_integer(dspfvs) and _positive_integer(decim):
        delay = bruker_dsp_table.get(dspfvs, {}).get(decim)
        if delay is not None:
            return DelayDecision(delay, "bruker_dsp_table", source)
    return DelayDecision(None, "unresolved", source,
                         "No positive finite GRPDLY or bruker_dsp_table entry for DSPFVS/DECIM")


def _jeol_delay(metadata):
    flag = metadata.get("digital_filter")
    if isinstance(flag, str):
        flag = {"TRUE": True, "FALSE": False}.get(flag)
    source = _source_fields(metadata, "digital_filter")
    if flag is False:
        return DelayDecision(0.0, "filter_disabled", source)
    if flag is not True:
        return DelayDecision(None, "unresolved", source, "digital_filter must be TRUE or FALSE")

    keys = ("digital_filter", "orders", "factors")
    if "filter_factor" in metadata:
        keys += ("filter_factor",)
    source = _source_fields(metadata, *keys)
    orders, factors = metadata.get("orders"), metadata.get("factors")
    if not isinstance(orders, (tuple, list)) or not orders or not all(map(_positive_integer, orders)):
        return DelayDecision(None, "unresolved", source, "orders must contain positive integers")
    if orders[0] != len(orders) - 1:
        return DelayDecision(None, "unresolved", source, "orders[0] must equal the number of stage orders")
    if not isinstance(factors, (tuple, list)) or not all(map(_positive_integer, factors)):
        return DelayDecision(None, "unresolved", source, "factors must contain positive integers")
    if len(factors) != orders[0]:
        return DelayDecision(None, "unresolved", source, "orders and factors must describe the same number of stages")

    decimation = prod(factors)
    filter_factor = metadata.get("filter_factor")
    if "filter_factor" in metadata and (not _positive_integer(filter_factor) or filter_factor != decimation):
        return DelayDecision(None, "unresolved", source, "filter_factor must equal the product of factors")

    # jeolconverter and OpenReadout express each stage delay in output dwell times.
    delay = 0.0
    for order, factor in zip(orders[1:], factors):
        delay += (order - 1) / decimation
        decimation //= factor
    return DelayDecision(0.5 * delay, "jeol_fir_cascade", source)


def resolve_delay(raw: RawFID) -> DelayDecision:
    """Resolve the delay from the vendor's digital filter metadata."""
    metadata = raw.digital_filter or {}
    if raw.vendor == "bruker":
        return _bruker_delay(metadata)
    if raw.vendor == "jeol":
        return _jeol_delay(metadata)
    if raw.vendor == "synthetic" and raw.digital_filter is None:
        return DelayDecision(0.0, "synthetic_no_filter", (("digital_filter", None),))
    return DelayDecision(None, "unresolved", (), f"Unsupported vendor {raw.vendor!r}")


def remove_delay(raw: RawFID, decision: DelayDecision) -> DelayRemoval:
    """Advance a FID by the resolved delay and trim the wrapped tail."""
    delay = decision.delay_points
    if delay is None:
        raise ValueError(f"Cannot remove unresolved digital filter delay: {decision.reason}")
    if not isfinite(delay) or delay < 0:
        raise ValueError("Digital filter delay must be finite and nonnegative")
    if delay == 0:
        return DelayRemoval(raw.fid.copy(), decision, "bypass", 0, "unchanged")

    # Keep only samples whose shifted time lies within the acquisition.
    tail_trim = ceil(delay)
    if tail_trim >= len(raw.fid):
        raise ValueError("Digital filter delay leaves no usable FID points")

    # Signed frequencies give one continuous phase ramp across the carrier.
    frequencies = np.fft.fftfreq(len(raw.fid))
    advance = np.exp(2j * np.pi * frequencies * delay)
    shifted = np.fft.ifft(np.fft.fft(raw.fid) * advance)
    return DelayRemoval(shifted[:-tail_trim], decision, "numpy.fft.fft / numpy.fft.ifft",
                        tail_trim, "unweighted after shift")
