"""Read time-domain Bruker acquisitions."""

import re
from math import isfinite
from os import fstat
from pathlib import Path

from nmrglue.fileio import bruker

from .model import UnsupportedInput, _raw_fid, _require_proton, _validate_fields

FILTER_KEYS = ("GRPDLY", "DECIM", "DSPFVS")
JCAMP_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[Ee][+-]?[0-9]+)?")


class _ParameterLines:
    # nmrglue 0.12 waits for the rest of a string or array even at EOF.
    def __init__(self, lines):
        self.lines = lines
        self.next_line = next(lines, None)

    def __iter__(self):
        return self

    def __next__(self):
        if self.next_line is None:
            raise StopIteration
        line = self.next_line
        self.next_line = next(self.lines, None)
        return line

    def readline(self):
        if self.next_line is None or self.next_line.startswith("##"):
            raise EOFError("incomplete JCAMP parameter")
        return next(self)


def _parameter_one(line, lines):
    declaration = line.partition("=")[2].partition("$$")[0]
    bounds = re.fullmatch(r"\s*\(\s*(\d+)\s*\.\.\s*(\d+)\s*\)(.*)", declaration)
    if bounds is None:
        raise ValueError("invalid JCAMP array range")
    first, last = int(bounds[1]), int(bounds[2])
    tokens = bounds[3].split()
    # nmrglue stops once it has enough tokens, leaving later surplus lines unread.
    while lines.next_line is not None and not lines.next_line.startswith("##"):
        tokens.extend(next(lines).partition("$$")[0].split())
    if last < first or len(tokens) != last - first + 1:
        raise ValueError("JCAMP array length does not match its range")
    if not all(JCAMP_NUMBER.fullmatch(token) for token in tokens):
        raise ValueError("invalid JCAMP array number")
    values = [bruker.parse_jcamp_value(token.lower()) for token in tokens]
    if not all(type(value) in (int, float) and isfinite(value) for value in values):
        raise ValueError("JCAMP array values must be finite numbers")
    if first <= 1 <= last:
        return values[1 - first]
    return None


def _decode_acqus(buffer):
    buffer = buffer.removeprefix(b"\xef\xbb\xbf")
    errors = []
    for encoding in ("utf-8", "cp1252"):
        try:
            return buffer.decode(encoding)
        except UnicodeDecodeError as error:
            errors.append(error)
    # Diagnose the error after the longest prefix either decoder accepted.
    raise max(errors, key=lambda error: error.start)


def _undecodable_parameter(error):
    encoding = "utf-8" if error.encoding == "utf-8" else "cp1252"
    prefix = error.object[:error.start].decode(encoding)
    for line in reversed(prefix.splitlines(keepends=True)):
        if line.startswith("##$"):
            key, separator, _ = line[3:].partition("=")
            return key.strip() if separator else None
        if line.startswith("##") or (line.lstrip().startswith("$$") and not line.endswith(("\n", "\r"))):
            break
    return None


def _read_acqus(path):
    try:
        buffer = path.read_bytes()
    except FileNotFoundError as error:
        raise UnsupportedInput("missing_acqus", path.parent, "acqus file not found") from error
    key = None
    try:
        text = _decode_acqus(buffer)
        lines = _ParameterLines(line.rstrip() for line in text.splitlines())
        params = {}
        for line in lines:
            if line.startswith("##END="):
                return params
            if line.startswith("##$"):
                key = line[3:].partition("=")[0].strip()
                if key in ("D", "P"):
                    params[key] = _parameter_one(line, lines)
                else:
                    if key == "NUC1":
                        line = line.partition("$$")[0].rstrip()
                    key, value = bruker.parse_jcamp_line(line, lines)
                    params[key] = value
    except (EOFError, ValueError, OverflowError) as error:
        if isinstance(error, UnicodeDecodeError):
            key = _undecodable_parameter(error)
        detail = "invalid or incomplete JCAMP parameter" + (f" {key}" if key else "")
        raise UnsupportedInput("truncated", path, detail) from error
    raise UnsupportedInput("truncated", path, "missing JCAMP ##END= line")


def _read_fid(path, params):
    td = params.get("TD")
    byte_order = params.get("BYTORDA")
    data_type = params.get("DTYPA")
    if type(td) is not int or td <= 0 or td % 2:
        raise UnsupportedInput("truncated", path, "TD must give a positive, even scalar count")
    for key, value, allowed in (("BYTORDA", byte_order, (0, 1)), ("DTYPA", data_type, (0, 2))):
        if value not in allowed:
            raise UnsupportedInput("truncated", path, f"missing or invalid {key}")
    is_float = data_type == 2
    scalar_bytes = 8 if is_float else 4
    expected = td * scalar_bytes
    padded = ((expected + 1023) // 1024) * 1024
    try:
        stream = path.open("rb")
    except FileNotFoundError as error:
        raise UnsupportedInput("missing_fid", path.parent, "fid file not found") from error
    with stream:
        size = fstat(stream.fileno()).st_size
        if size > padded:
            raise UnsupportedInput("arrayed_or_2d", path, "fid exceeds the size of a single trace")
        if size not in (expected, padded):
            raise UnsupportedInput("truncated", path, f"expected {expected} fid bytes, found {size}")
        try:
            data = bruker.get_data(stream, big=byte_order == 1, isfloat=is_float)
        except ValueError as error:
            raise UnsupportedInput("truncated", path, "incomplete binary scalar") from error
    if data.nbytes != size:
        raise UnsupportedInput("truncated", path, "fid byte size changed while reading")
    return bruker.complexify_data(data[:td])


def _dimension_evidence(folder, params):
    mode, blocks = params.get("PARMODE"), params.get("NBL")
    if mode not in (None, 0):
        return f"PARMODE={mode}"
    if params.get("FnTYPE") == 2:
        return "FnTYPE=2"
    if blocks not in (None, 0, 1):
        return f"NBL={blocks}"
    return next((name for name in ("acqu2s", "acqu3s", "acqu4s")
                 if (folder / name).is_file()), None)


def read_bruker(folder):
    """Read one complex proton fid using its acquisition status file."""
    folder = Path(folder)
    if not folder.is_dir():
        raise UnsupportedInput("not_found", folder, "Bruker acquisition directory not found")
    acqus = folder / "acqus"
    fid_path = folder / "fid"
    if not acqus.is_file():
        raise UnsupportedInput("missing_acqus", folder, "acqus file not found")
    has_fid = fid_path.is_file()
    has_ser = (folder / "ser").is_file()
    if not (has_fid or has_ser):
        raise UnsupportedInput("missing_fid", folder, "fid file not found")
    params = _read_acqus(acqus)
    if "SW_h" in params and params["SW_h"] is None:
        raise UnsupportedInput("truncated", acqus, "empty SW_h")
    for key in ("D", "P"):
        params[f"{key}[1]"] = params.get(key)
    _validate_fields(
        acqus, params, integers=("TD", "NS", "DS", "BYTORDA", "DTYPA", "AQ_mod", "PARMODE",
                                 "FnTYPE", "NBL", "NC", "DECIM", "DSPFVS"),
        reals=("SW_h", "SFO1", "BF1", "O1", "D[1]", "P[1]", "RG", "TE", "GRPDLY"),
        strings=("NUC1", "PULPROG"),
    )
    evidence = _dimension_evidence(folder, params)
    if evidence or (has_fid and has_ser):
        detail = (f"acqus indicates multiple traces or dimensions ({evidence})" if evidence
                  else "conflicting acquisition files")
        if has_ser:
            detail = f"found {'fid and ser' if has_fid else 'ser'}; {detail}"
        raise UnsupportedInput("arrayed_or_2d", folder, detail)
    if has_ser:
        raise UnsupportedInput("ser_not_fid", folder, "found ser without fid or dimension evidence")
    mode = params.get("AQ_mod")
    if mode not in (0, 1, 2, 3):
        raise UnsupportedInput("truncated", folder, "missing or invalid AQ_mod")
    if mode in (0, 2):
        raise UnsupportedInput("arrayed_or_2d", folder, f"AQ_mod={mode} is not complex")
    nucleus = params.get("NUC1") or ""
    if nucleus.strip().casefold() == "off":
        nucleus = ""
    _require_proton(folder, "NUC1", nucleus, ("1h", "¹h"))

    offset = params.get("O1")
    base = params.get("BF1")
    carrier = None if offset is None or not base else offset / base
    if carrier is not None and not isfinite(carrier):
        raise UnsupportedInput("truncated", acqus, "carrier ppm is out of range")
    fid = _read_fid(fid_path, params)
    digital_filter = {key: params[key] for key in FILTER_KEYS if key in params}
    nc = params.get("NC")
    scale = "NC is absent" if nc is None else f"NC={nc}"
    return _raw_fid(
        "bruker", folder, fid, f"Stored Bruker values ({scale}). No NC correction or scan normalization applied.",
        {
            "nucleus": ("1H", "NUC1"),
            "observe_frequency_mhz": (params.get("SFO1"), "SFO1"),
            "spectral_width_hz": (params.get("SW_h"), "SW_h"),
            "carrier_ppm": (carrier, "O1 / BF1"),
            "complex_points": (len(fid), "TD"),
            "relaxation_delay_s": (params["D[1]"], "D[1]"),
            "pulse_width_us": (params["P[1]"], "P[1]"),
            "flip_angle_deg": (None, None),
            "pulse_90_us": (None, None),
            "scans": (params.get("NS"), "NS"),
            "dummy_scans": (params.get("DS"), "DS"),
            "receiver_gain": (params.get("RG"), "RG"),
            "temperature_k": (params.get("TE"), "TE"),
            "pulse_program": (params.get("PULPROG"), "PULPROG"),
            "digital_filter": (digital_filter or None, ", ".join(digital_filter)),
        },
        intensity_scale=1.0,
        storage_scale=None if nc is None else scale,
    )
