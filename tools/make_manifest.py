"""Verify local qNMR archives and inventory their extracted files."""

import argparse
import csv
import hashlib
import json
import os
import re
import stat
import sys
import tempfile
from collections import Counter
from pathlib import Path

COLUMNS = ("deposit_id", "archive", "path", "size_bytes", "sha256", "kind", "licence", "doi")
BRUKER_NAMES = {"fid", "ser", "acqus", "acqu2s", "procs", "1r"}
WORKBOOK_SUFFIXES = {".xlsx", ".xls", ".numbers", ".csv"}


def require_text(value, field):
    if not isinstance(value, str) or not value.strip() or "\0" in value:
        raise ValueError(f"{field} must be a non-empty string without NULs")
    try:
        value.encode("utf-8")
    except UnicodeEncodeError as error:
        raise ValueError(f"{field} must be valid UTF-8") from error


def require_path_component(value, field):
    require_text(value, field)
    if value == ".." or Path(value).name != value:
        raise ValueError(f"{field} must be a single file or directory name")


def validate_deposits(deposits):
    if not isinstance(deposits, list) or not deposits:
        raise ValueError("deposits must be a non-empty list")
    for index, deposit in enumerate(deposits, 1):
        label = f"deposit {index}"
        if not isinstance(deposit, dict):
            raise ValueError(f"{label} must be an object")
        require_path_component(deposit.get("id"), f"{label} id")
        for field in ("doi", "licence"):
            require_text(deposit.get(field), f"{label} {field}")
        archives = deposit.get("archives")
        if not isinstance(archives, list) or not archives:
            raise ValueError(f"{label} archives must be a non-empty list")
        for number, archive in enumerate(archives, 1):
            archive_label = f"{label} archive {number}"
            if not isinstance(archive, dict):
                raise ValueError(f"{archive_label} must be an object")
            require_path_component(archive.get("local_name"), f"{archive_label} local_name")
            md5 = archive.get("md5")
            if not isinstance(md5, str) or not re.fullmatch(r"[0-9a-fA-F]{32}", md5):
                raise ValueError(f"{archive_label} md5 must contain 32 hexadecimal digits")
            size = archive.get("size_bytes")
            if type(size) is not int or size < 0:
                raise ValueError(f"{archive_label} size_bytes must be a non-negative integer")


def extracted_folder(archive_path):
    return archive_path.with_name(f"x_{archive_path.stem}")


def validate_output(out, deposits_file, archive_paths):
    resolved = out.resolve()
    if resolved == deposits_file.resolve():
        raise ValueError("--out matches the --deposits file")
    destination = out.parent.resolve() / out.name
    for archive_path in archive_paths:
        if resolved == archive_path.resolve():
            raise ValueError(f"--out matches archive {archive_path}")
        folder = extracted_folder(archive_path).resolve()
        if resolved.is_relative_to(folder) or destination.is_relative_to(folder):
            raise ValueError(f"--out is inside extracted folder {folder}")


def checksum(path, digest):
    size = 0
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def check_archive(path, archive):
    label = f"{path.parent.name}/{path.name}"
    try:
        if not stat.S_ISREG(path.lstat().st_mode):
            print(f"MISMATCH {label} (expected a regular file)")
            return False
        size, md5 = checksum(path, hashlib.md5(usedforsecurity=False))
    except FileNotFoundError:
        print(f"MISSING {label}")
        return False
    except OSError as error:
        print(f"MISMATCH {label} ({error})")
        return False

    if size != archive["size_bytes"] or md5 != archive["md5"].lower():
        print(f"MISMATCH {label} (size {size}, MD5 {md5})")
        return False
    print(f"OK {label}")
    return True


def file_kind(path):
    if path.name == ".DS_Store" or path.name.startswith("._") or "__MACOSX" in path.parts[:-1]:
        return "os_junk"
    if path.name in BRUKER_NAMES:
        return f"bruker_{path.name}"
    if path.suffix.lower() == ".jdf":
        return "jdf"
    if path.suffix.lower() in WORKBOOK_SUFFIXES:
        return "workbook"
    return "other"


def raise_walk_error(error):
    raise error


def extracted_rows(folder, deposit, archive_name):
    if not stat.S_ISDIR(folder.lstat().st_mode):
        raise ValueError(f"expected an extracted directory at {folder}")
    for directory, _, filenames in os.walk(folder, onerror=raise_walk_error):
        for filename in filenames:
            path = Path(directory) / filename
            if stat.S_ISREG(path.lstat().st_mode):
                relative = path.relative_to(folder)
                size, sha256 = checksum(path, hashlib.sha256())
                yield (deposit["id"], archive_name, relative.as_posix(), size, sha256,
                       file_kind(relative), deposit["licence"], deposit["doi"])


def write_manifest(out, rows):
    stream = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=out.parent, delete=False
    )
    temporary = Path(stream.name)
    try:
        with stream:
            csv.writer(stream, lineterminator="\n").writerow(COLUMNS)
            # Quote paths with embedded carriage returns as well as commas and newlines.
            csv.writer(stream, lineterminator="\n", quoting=csv.QUOTE_ALL).writerows(rows)
        temporary.replace(out)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="directory containing deposit IDs")
    parser.add_argument("--deposits", required=True, type=Path, help="deposit metadata JSON")
    parser.add_argument("--out", required=True, type=Path, help="output CSV")
    args = parser.parse_args()

    try:
        deposits = json.loads(args.deposits.read_text(encoding="utf-8"))
        validate_deposits(deposits)
        archives = [
            (deposit, archive, args.root / deposit["id"] / archive["local_name"])
            for deposit in deposits
            for archive in deposit["archives"]
        ]
        validate_output(args.out, args.deposits, (path for _, _, path in archives))
        checks = [check_archive(path, archive) for _, archive, path in archives]
        if not all(checks):
            return 1

        rows = []
        for deposit, _, path in archives:
            rows.extend(extracted_rows(extracted_folder(path), deposit, path.name))
        rows.sort(key=lambda row: (row[0], row[2], row[1]))
        write_manifest(args.out, rows)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"make_manifest: {error}", file=sys.stderr)
        return 1

    counts = {deposit["id"]: Counter() for deposit in deposits}
    for deposit_id, _, _, _, _, kind, _, _ in rows:
        counts[deposit_id][kind] += 1
    for deposit_id, kinds in sorted(counts.items()):
        details = ", ".join(f"{kind}={count}" for kind, count in sorted(kinds.items()))
        print(f"{deposit_id}: {sum(kinds.values())} files {details}".rstrip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
