"""Build the development case manifest from inventory and acquisition metadata."""

import argparse
import csv
import hashlib
import json
import re
import struct
import sys
import tempfile
from os import fstat
from pathlib import Path, PurePosixPath

PROJECT = Path(__file__).resolve().parents[1]
NISHIZAKI = "nishizaki_QFY6L9"
ARCHIVE_STEM = "334_Dataverse_20200705"
STORED_ROOT = ".cache/veri/"
FID02 = "FID02_Figure2/01_5mmNMRtube/"
FID05 = "FID05_Figure4/01_13C_decoupled_Aquisition/N1/Delay60sec/"
FID06 = "FID06_Table4/01_Isoniazid/"
PILOT_PATHS = {
    "FID02": (
        FID02 + "Pattern4/0724_2019/04_MAandCAF_T1A_5mm_6435u15130u9000uL_yuzo190724_RG30_D20s_qH-1-1.jdf",
        FID02 + "Pattern4/0724_2019/01_MAandDMSO2_T1A_5mm_4226u1534u6000uL_yuzo190724_RG30_D20s_qH-1-1.jdf",
    ),
    "FID05": (
        FID05 + "13_DMSO2_CD3OD_5mm_1728u2000uL_yuzo190815_qH_D60s_NS16-1-1.jdf",
        FID05 + "12_DMSO2_DMSO_5mm_1599u1500uL_yuzo190815_qH_D60s_NS16-1-1.jdf",
    ),
    "FID06": (
        FID06 + "BTMSBandINH_5mm_1064u8463u1500uL_yuzo190817_qH_D60s_NS16-1-1.jdf",
        FID06 + "13_DMSO2_CD3OD_5mm_1728u2000uL_yuzo190817_qH_D60s_NS16-1-1.jdf",
    ),
}
DATES_01 = {
    "Pattern3": {"0711_2019", "0725_2019", "0727_2019"},
    "Pattern4": {"0711_2019", "0718_2019", "0719_2019", "0724_2019"},
}
COHORT_HASHES = {
    "01": "9bfca77ca06a16ab62026940962503a03758b0d1f3288b773fd2059fcd0f8574",
    "09": "19503cff5ab8be126deeb5751dc97e6b3eb22523c7a97fc1d198c3e4675f0510",
    "July": "11f6ef9c493117697469de9ac2b1c030683de1a9a36f064ee6ffd65bea86e8a3",
}
WINDOW_FIELDS = ("path", "bounds", "digests", "settings")
INTEGRAL_CASES = "tests/fixtures/integral_cases.json"
EXPOSURES = [
    {"source": "job 07", "scope": "FID03 Sheet1 column E area cells for solutions 09–11"},
    {"source": "job 07", "scope": "FID02 5 mm sheet rows 146–149, including L147 and L148"},
]


def _archive_key(deposit, archive):
    return f"{deposit}/x_{PurePosixPath(archive).stem}"


def _local_key(key):
    deposit, archive, path = key
    return f"{_archive_key(deposit, archive)}/{path}"


def _nishizaki_key(path):
    return NISHIZAKI, ARCHIVE_STEM + ".zip", ARCHIVE_STEM + "/" + path


FID02_PREFIX = STORED_ROOT + _local_key(_nishizaki_key(FID02))
CONTROLS_09 = {
    "tests/fixtures/jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf",
    "tests/fixtures/bruker/zenodo17083860_833-1/10",
    *(STORED_ROOT + _local_key(_nishizaki_key(path))
      for figure in ("FID05", "FID06") for path in PILOT_PATHS[figure]),
}


def _validate_key(key):
    deposit, archive, path = key
    if (not all(isinstance(part, str) for part in key) or "/" in deposit or "/" in archive
            or any(part in {"", ".", ".."} or "\\" in part or ":" in part or "\x00" in part
                   for part in (deposit, archive, *path.split("/")))):
        raise ValueError(f"Invalid inventory key {key!r}")


def _metadata_path(data_root, key, is_jdf):
    deposit, archive, path = key
    outside_archive = f"Metadata path leaves archive folder for inventory key {key!r}"
    try:
        root = data_root.resolve()
        folder = (root / _archive_key(deposit, archive)).resolve()
        if folder == root or not folder.is_relative_to(root):
            raise ValueError(outside_archive)
        metadata = folder / path
        if not is_jdf:
            metadata = metadata.with_name("acqus")
        metadata = metadata.resolve()
    except (OSError, RuntimeError) as error:
        raise ValueError(f"Cannot resolve metadata path for inventory key {key!r}") from error
    if metadata == folder or not metadata.is_relative_to(folder):
        raise ValueError(outside_archive)
    return metadata


def _acquisition_id(path, is_jdf):
    try:
        if not is_jdf:
            with path.open("rb") as stream:
                for line in stream:
                    match = re.match(rb"##\$DATE\s*=\s*([+-]?\d+)\s*(?:\$\$.*)?$", line)
                    if match:
                        return {"date": int(match[1])}
            return None
        with path.open("rb") as stream:
            from nmrglue.fileio import jeol
            from spinivo.readers.jeol import _ConversionTable, _read_parameters

            # Read the fixed JDF header and the parameter block only.
            buffer, header = jeol.read_header(jeol.IOBuffer(stream.read(1360), _ConversionTable))
            start, length = header["param_start"], header["param_length"]
            if (header["file_identifier"] != "JEOL.NMR" or start < buffer.position
                    or start + length > min(header["data_start"], fstat(stream.fileno()).st_size)):
                raise ValueError(f"Invalid JDF parameter block in {path}")
            stream.seek(start)
            params = jeol.IOBuffer(stream.read(length), _ConversionTable)
            values = _read_parameters(params, 0, length, header["endian"], path)
            if {"actual_start_time", "end_time"} <= values.keys():
                return {"start": values["actual_start_time"], "end": values["end_time"], "title": header["title"]}
    except FileNotFoundError:
        pass
    except (struct.error, UnicodeError, KeyError, TypeError, OverflowError) as error:
        raise ValueError(f"Invalid JDF parameter block in {path}") from error
    return None


def _preparation(name):
    recipe = re.search(r"_((?:\d+u)+\d+uL)_", name)
    if recipe is None:
        return None
    label = name.split("_", 1)[0]
    return {"id": f"{label}:{recipe[1]}", "label": label, "recipe_token": recipe[1], "source": "Nishizaki filename"}


def _links(preparations):
    shared = "FID02/FID04 filenames: same masses, one third volume"
    pairs = [
        ("shared_weighing", "01:4226u1534u6000uL", "05:4226u1534u2000uL", shared),
        ("shared_weighing", "02:3913u2894u6000uL", "06:3913u2894u2000uL", shared),
        ("shared_weighing", "03:4280u5692u6000uL", "07:4280u5692u2000uL", shared),
        ("shared_weighing", "04:6435u15130u9000uL", "08:6435u15130u3000uL", shared),
        ("preparation_alias", "MAandCAF:6435u15130u9000uL", "04:6435u15130u9000uL",
         "FID01 MAandCAF and FID02 solution 04 filenames: same recipe"),
        ("preparation_alias", "12:2346u2000uL", "16:2346u2000uL",
         "FID05 solution 12 D2O nutation and solution 16 filenames: same recipe and solvent"),
    ]
    for _, left, right, _ in pairs:
        for identity in (left, right):
            if identity not in preparations:
                raise ValueError(f"Missing preparation {identity} for link {left} / {right}")
    return [
        {"kind": kind, "preparation_ids": [left, right], "source": source, "status": "possible"}
        for kind, left, right, source in pairs
    ]


def _runs(rows, data_root):
    runs, preparations, by_hash = {}, {}, {}
    for key in rows:
        _validate_key(key)
    for key, row in sorted(rows.items()):
        deposit, archive, path = key
        is_jdf = row["kind"] == "jdf"
        if is_jdf:
            stem = path[:-8]
            if path.lower().endswith("-1-2.jdf") and (deposit, archive, stem + "-1-1" + path[-4:]) in rows:
                continue
            processed = stem + "-1-2" + path[-4:] if path.lower().endswith("-1-1.jdf") else None
        elif row["kind"] in {"bruker_fid", "bruker_ser"}:
            folder = PurePosixPath(path).parent
            if row["kind"] == "bruker_ser" and (deposit, archive, str(folder / "fid")) in rows:
                continue
            processed = str(folder / "pdata/1/1r")
        else:
            continue
        processed_key = deposit, archive, processed
        prep = _preparation(PurePosixPath(path).name) if is_jdf and deposit == NISHIZAKI else None
        if prep:
            preparations[prep["id"]] = prep
        run_id = "/".join(key)
        digest = row["sha256"]
        runs[key] = {
            "id": run_id, "deposit": deposit, "raw_key": list(key),
            "processed_key": list(processed_key) if processed_key in rows else None,
            "sha256": digest, "acquisition_id": _acquisition_id(_metadata_path(data_root, key, is_jdf), is_jdf),
            "alias_of": by_hash.get(digest, {}).get("id"),
            "preparation_id": prep["id"] if prep else None, "fixture_aliases": [],
        }
        by_hash.setdefault(digest, runs[key])
    return runs, preparations, by_hash


def _cohort_names(path):
    parts = path.removeprefix(FID02_PREFIX).split("/") if path.startswith(FID02_PREFIX) else []
    names = []
    if len(parts) == 3 and parts[2].split("_", 1)[0] in {"01", "02", "03", "04"}:
        pattern, date, name = parts
        if date in DATES_01.get(pattern, ()) and name.endswith("D20s_qH-1-1.jdf"):
            names.append("01")
        if pattern == "Pattern4" and name.endswith("_qH-1-1.jdf"):
            names.append("09")
            if date.startswith("07"):
                names.append("July")
    if path in CONTROLS_09:
        names.append("09")
    return names


def _window_cases(windows):
    cases = []
    members = {name: [] for name in COHORT_HASHES}
    for window in windows:
        fields = {name: window[name] for name in WINDOW_FIELDS}
        if not isinstance(fields["path"], str):
            raise ValueError("Baseline window path must be a string")
        identity = json.dumps(list(fields.values()), sort_keys=True, ensure_ascii=True, separators=(", ", ": "))
        case_id = "window-" + hashlib.sha256(identity.encode("utf-8")).hexdigest()
        names = _cohort_names(fields["path"])
        for name in names:
            members[name].append((identity, case_id))
        cases.append(dict(fields, id=case_id, kind="window", cohorts=names, split="development"))
    cohorts = {}
    for name, entries in members.items():
        ordered = sorted(entries)
        digest = hashlib.sha256("\n".join(identity for identity, _ in ordered).encode("utf-8")).hexdigest()
        if digest != COHORT_HASHES[name]:
            raise ValueError(f"Cohort {name} hash mismatch: expected {COHORT_HASHES[name]}, got {digest}")
        cohorts[name] = {"members": [case_id for _, case_id in ordered], "sha256": digest}
    return cases, cohorts


def _attach_windows(cases, runs, by_hash, rows):
    by_path = {STORED_ROOT + _local_key(key): run for key, run in runs.items()}
    by_bruker_digests = {}
    for (deposit, archive, path), run in runs.items():
        if rows[deposit, archive, path]["kind"] == "bruker_fid":
            acqus = rows.get((deposit, archive, str(PurePosixPath(path).with_name("acqus"))))
            if acqus:
                by_bruker_digests.setdefault((run["sha256"], acqus["sha256"]), run)
    for case in cases:
        path = case["path"]
        name = PurePosixPath(path).name
        if path.startswith("tests/fixtures/"):
            if name.endswith(".jdf"):
                digest = case["digests"][name]
                run = by_hash.get(digest)
            else:
                digest = case["digests"]["fid"]
                run = by_bruker_digests.get((digest, case["digests"]["acqus"]))
            if run and path not in run["fixture_aliases"]:
                run["fixture_aliases"].append(path)
        else:
            digest = case["digests"][name]
            run = by_path.get(path)
        if run is None or run["sha256"] != digest:
            raise ValueError(f"No matching inventory run for window {path}")
        case["run_ids"] = [run["id"]]


def build(manifest_csv, windows_json, data_root):
    """Return metadata for the frozen baseline, overlap case and pilot pairs."""
    windows = json.loads(Path(windows_json).read_text(encoding="utf-8"))
    cases, cohorts = _window_cases(windows)
    with Path(manifest_csv).open(encoding="utf-8", newline="") as stream:
        rows = {(row["deposit_id"], row["archive"], row["path"]): row for row in csv.DictReader(stream)}
    runs, preparations, by_hash = _runs(rows, Path(data_root))
    _attach_windows(cases, runs, by_hash, rows)

    integral_cases = json.loads((PROJECT / INTEGRAL_CASES).read_text(encoding="utf-8"))["cases"]
    overlap = next(case for case in integral_cases if case["id"] == "fid05_widened_dmso2")
    overlap_key = NISHIZAKI, ARCHIVE_STEM + ".zip", overlap["source"]
    cases.append({"id": overlap["id"], "kind": "overlap", "source": INTEGRAL_CASES,
                  "bounds": overlap["bounds"], "run_ids": [runs[overlap_key]["id"]],
                  "split": "development", "cohorts": []})
    pilot_pairs = []
    for figure, paths in PILOT_PATHS.items():
        pair = {"id": f"pilot-{figure}", "source": "docs/pilot-pairs.md"}
        for role, path in zip(("sample", "ec"), paths):
            run = runs[_nishizaki_key(path)]
            pair[role] = {"run_id": run["id"], "raw_key": run["raw_key"]}
        pilot_pairs.append(pair)
        cases.append({"id": pair["id"], "kind": "pilot_pair", "split": "development", "cohorts": [],
                      "run_ids": [pair["sample"]["run_id"], pair["ec"]["run_id"]]})
    return {"schema_version": 1, "runs": list(runs.values()),
            "preparations": [preparations[key] for key in sorted(preparations)],
            "links": _links(preparations), "cases": cases, "exposures": [row.copy() for row in EXPOSURES],
            "pilot_pairs": pilot_pairs, "cohorts": cohorts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=PROJECT / "data/manifest.csv")
    parser.add_argument("--windows", type=Path, default=PROJECT / "data/baseline_windows.json")
    parser.add_argument("--data-root", type=Path, default=Path(".cache/veri"))
    parser.add_argument("--out", type=Path, default=PROJECT / "data/case_manifest.json")
    args = parser.parse_args()
    resolved = args.out.resolve()
    if resolved in {args.manifest.resolve(), args.windows.resolve(), PROJECT / INTEGRAL_CASES}:
        parser.error("--out must differ from the input files")
    root = args.data_root.resolve()
    if resolved.is_relative_to(root) or (args.out.parent.resolve() / args.out.name).is_relative_to(root):
        parser.error("--out must be outside --data-root")
    result = build(args.manifest, args.windows, args.data_root)
    stream = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=args.out.parent, delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
        temporary.replace(args.out)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.path.insert(0, str(PROJECT / "src"))
    main()
