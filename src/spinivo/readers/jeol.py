"""Read time-domain JEOL JDF acquisitions."""

import struct
from math import degrees
from pathlib import Path

from nmrglue.fileio import jeol

from .model import UnsupportedInput, _raw_fid, _require_proton, _validate_fields

FILTER_KEYS = ("digital_filter", "orders", "factors", "filter_factor")
UNITS = {
    "relaxation_delay": ("second", 0), "x_pulse": ("second", 6), "x_90_width": ("second", 6),
    "x_sweep": ("hertz", 0), "x_freq": ("hertz", -6), "temp_get": ("kelvin", 0), "x_angle": ("degree", 0),
    "x_offset": ("ppm", 0),
}


class _ConversionTable(jeol.ConversionTable):
    # JDF encodes signed prefix codes in a four-bit field.
    prefix = {code & 15: name for code, name in jeol.ConversionTable.prefix.items() if name != "None"}


def _read_description(buffer, path):
    stream, header = jeol.read_header(jeol.IOBuffer(buffer, _ConversionTable))
    start, length = header["param_start"], header["param_length"]
    if start < stream.position or length < 16 or start + length > len(buffer):
        raise UnsupportedInput("truncated", path, "JDF parameter section is out of bounds")
    order = "<" if header["endian"] == "little_endian" else ">"
    size, _, count, total = struct.unpack_from(order + "4I", buffer, start)
    if size != 64 or not count * size <= total <= length - 16:
        raise UnsupportedInput("truncated", path, "JDF parameter records are out of bounds")
    stages = {}
    # nmrglue multiplies even strings by 10 ** unit_scaler.
    for record in range(start + 16, start + 16 + count * size, size):
        scaler = struct.unpack_from(order + "h", buffer, record + 4)[0]
        value_type = struct.unpack_from(order + "i", buffer, record + 32)[0]
        if abs(scaler) > 308 or (value_type == 0 and scaler):
            raise UnsupportedInput("truncated", path, "invalid JDF parameter scaler")
        name = buffer[record + 36:record + size].replace(b"\x00", b"").replace(b" ", b"").lower()
        if name in (b"orders", b"factors"):
            if value_type != 0:
                raise UnsupportedInput("truncated", path, f"expected a JDF stage string for {name.decode()}")
            # Preserve stage boundaries removed by nmrglue's string parser.
            tokens = buffer[record + 16:record + 32].split(b"\x00", 1)[0].split()
            try:
                stages[name.decode()] = tuple(map(int, tokens))
            except ValueError as error:
                raise UnsupportedInput("truncated", path, f"invalid JDF stages for {name.decode()}") from error
    _, params = jeol.read_parameters(stream, start, header["endian"])
    params.update(stages)
    return {"header": header, "parameters": params}


def _convert_units(params, path):
    offset_source = "x_offset"
    for key, (expected, target_exponent) in UNITS.items():
        if key not in params:
            continue
        (prefix, power, unit), *extra = params["info"][key]["units"]
        if power != 1 or any(p or base for _, p, base in extra):
            raise UnsupportedInput("truncated", path, f"invalid unit for {key}")
        exponent = jeol.ConversionTable.unit_prefix_table[prefix.title()]
        value = params[key] * 10.0 ** (exponent + target_exponent)
        if key == "temp_get" and unit == "celsius":
            value += 273.15
        elif key == "x_angle" and unit == "radian":
            value = degrees(value)
        elif key == "x_offset" and unit == "hertz":
            frequency = params.get("x_freq")  # Already converted to MHz.
            if frequency is not None and frequency <= 0:
                raise UnsupportedInput("truncated", path, "x_freq must be positive for an offset in Hz")
            value = None if frequency is None else value / frequency
            offset_source = "x_offset, x_freq"
        elif unit != expected:
            raise UnsupportedInput("truncated", path, f"invalid unit for {key}: {unit!r}")
        params[key] = value
    return offset_source


def read_jeol(path):
    """Read a 1D proton JDF and preserve the sample values returned by nmrglue."""
    path = Path(path)
    if not path.is_file():
        raise UnsupportedInput("not_found", path, "JDF file not found")
    try:
        buffer = path.read_bytes()
    except FileNotFoundError as error:
        raise UnsupportedInput("not_found", path, "JDF file not found") from error
    if not buffer.startswith(b"JEOL.NMR"):
        raise UnsupportedInput("truncated", path, "missing or incomplete JEOL file identifier")
    try:
        description = _read_description(buffer, path)
    except (struct.error, UnicodeError, KeyError, TypeError, OverflowError) as error:
        raise UnsupportedInput("truncated", path, "invalid or incomplete JDF header") from error

    header = description["header"]
    params = description["parameters"]
    if (header["data_format"] != "one_d" or header["data_dimension_number"] != 1
            or jeol.ndims(description) != 1 or any(n > 1 for n in header["data_points"][1:])):
        raise UnsupportedInput("arrayed_or_2d", path, f'JDF layout is {header["data_format"]}')
    axis = header["data_axis_type"][0]
    if axis == "real":
        raise UnsupportedInput("processed_spectrum", path, "JDF contains a real spectrum")
    if axis != "complex":
        raise UnsupportedInput("arrayed_or_2d", path, f"JDF axis {axis!r} is not complex")
    if header["units"][0][2] in ("ppm", "hertz"):
        raise UnsupportedInput("processed_spectrum", path, "JDF has a frequency-domain axis")
    integers = ("x_points", "scans", "x_prescans")
    reals = (*UNITS, "recvr_gain")
    strings = ("x_domain", "experiment")
    for key in (*integers, *reals, *strings):
        if key in params and params["info"][key]["value_type"] == "infinity":
            raise UnsupportedInput("truncated", path, f"invalid value for {key}")
    _validate_fields(path, params, integers=integers, reals=reals, strings=strings)
    offset_source = _convert_units(params, path)
    _validate_fields(path, params, reals=UNITS)
    _require_proton(path, "x_domain", params.get("x_domain"), ("proton", "1h", "h1", "¹h"))

    # JDF counts axis padding in stored_points and uses offsets to select samples.
    stored_points = header["data_points"][0]
    first = header["data_offset_start"][0]
    last = header["data_offset_stop"][0]
    point_count = last - first + 1
    scalar_bytes = {"float64": 8, "float32": 4}.get(header["data_type"])
    data_start = header["data_start"]
    data_length = header["data_length"]
    param_end = header["param_start"] + header["param_length"]
    if (scalar_bytes is None or not 0 <= first <= last < stored_points
            or data_start < param_end or data_length != 2 * stored_points * scalar_bytes
            or max(header["total_size"], data_start + data_length) > len(buffer)
            or params.get("x_points", point_count) != point_count):
        raise UnsupportedInput("truncated", path, "JDF point count or data bounds are inconsistent")
    data = jeol.read_bin_data(description, buffer)
    fid = jeol.truncate_data(description, jeol.reorganize_1d(description, data))
    if fid.shape != (point_count,):
        raise UnsupportedInput("truncated", path, f"expected {point_count} complex JDF samples")

    digital_filter = {key: params[key] for key in FILTER_KEYS if key in params}
    return _raw_fid(
        "jeol", path, fid, "Stored JEOL values. No additional intensity scaling or scan normalization applied.",
        {
            "nucleus": ("1H", "x_domain"),
            "observe_frequency_mhz": (params.get("x_freq"), "x_freq"),
            "spectral_width_hz": (params.get("x_sweep"), "x_sweep"),
            "carrier_ppm": (params.get("x_offset"), offset_source),
            "complex_points": (
                point_count, "x_points" if "x_points" in params else "data_offset_start, data_offset_stop",
            ),
            "relaxation_delay_s": (params.get("relaxation_delay"), "relaxation_delay"),
            "pulse_width_us": (params.get("x_pulse"), "x_pulse"),
            "flip_angle_deg": (params.get("x_angle"), "x_angle"),
            "pulse_90_us": (params.get("x_90_width"), "x_90_width"),
            "scans": (params.get("scans"), "scans"),
            "dummy_scans": (params.get("x_prescans"), "x_prescans"),
            "receiver_gain": (params.get("recvr_gain"), "recvr_gain"),
            "temperature_k": (params.get("temp_get"), "temp_get"),
            "pulse_program": (params.get("experiment"), "experiment"),
            "digital_filter": (digital_filter or None, ", ".join(digital_filter)),
        },
        intensity_scale=1.0,
        storage_scale="JEOL",
    )
