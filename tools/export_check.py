"""Dry run of a release export: check each file against data/redistribution.json.

From the repository, run python tools/export_check.py, or git ls-files -z | python tools/export_check.py -z.
Exit status 1 means an unknown, missing or changed file, or a cleared CC BY file without its credit files.
"""

import argparse
import fnmatch
import hashlib
import json
import os
import re
import sys
from pathlib import Path, PurePosixPath

PROJECT = Path(__file__).resolve().parents[1]
FOLDERS = ("src", "tests", "tools", "docs", "data")
TEXT_FIELDS = ("licence", "source", "modifications", "reason")
FAILURES = ("unknown", "missing", "mismatch", "attribution")


def load_allowlist(path):
    """Return the rules and the folder prefixes where only exact, hashed path rules apply."""
    allowlist = json.loads(Path(path).read_text(encoding="utf-8"))
    rules, exact_only = allowlist["rules"], tuple(allowlist.get("exact_only", ()))
    paths = set()
    for index, rule in enumerate(rules, 1):
        label = f"rule {index}"
        if ("path" in rule) == ("pattern" in rule):
            raise ValueError(f"{label} needs exactly one of path or pattern")
        for field in TEXT_FIELDS:
            if not isinstance(rule.get(field), str) or not rule[field].strip():
                raise ValueError(f"{label} {field} must be a non-empty string")
        if type(rule.get("distributable")) is not bool:
            raise ValueError(f"{label} distributable must be true or false")
        digest = rule.get("sha256")
        if digest is not None and ("pattern" in rule or not re.fullmatch(r"[0-9a-f]{64}", digest)):
            raise ValueError(f"{label} sha256 must be 64 lowercase hex digits on a path rule")
        if "path" in rule:
            if rule["path"] in paths:
                raise ValueError(f"{label} repeats path {rule['path']}")
            paths.add(rule["path"])
            if rule["path"].startswith(exact_only) and digest is None:
                raise ValueError(f"{label} needs sha256 under {exact_only}")
        credits = rule.get("credit_files", [])
        if "CC-BY-4.0" in rule["licence"] and (not credits or not all(isinstance(c, str) for c in credits)):
            raise ValueError(f"{label} CC-BY-4.0 needs credit_files")
    return rules, exact_only


def default_files(root):
    files = [path.name for path in root.iterdir() if path.is_file() and path.name != ".git"]
    for folder in FOLDERS:
        for path in (root / folder).rglob("*"):
            parts = path.relative_to(root).parts
            if path.is_file() and path.suffix != ".pyc" and not any(
                    part == "__pycache__" or part.endswith(".egg-info") for part in parts):
                files.append("/".join(parts))
    return sorted(files)


def read_names(stream, null):
    """Split NUL separated or line separated paths; line mode removes only the line ending."""
    data = stream.read()
    names = data.split(b"\0") if null else [line.removesuffix(b"\r") for line in data.split(b"\n")]
    return [os.fsdecode(name) for name in names if name]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def has_symlink(root, relative):
    path = root
    for part in relative.parts:
        path = path / part
        if path.is_symlink():
            return True
    return False


def check(root, names, rules, exact_only=()):
    """Return (status, path, licence, sha256, reason) per file; status is cleared, excluded or a failure."""
    exact = {rule["path"]: rule for rule in rules if "path" in rule}
    patterns = [rule for rule in rules if "pattern" in rule]
    results, credits = [], {}
    for name in names:
        relative = PurePosixPath(name)
        path = root / name
        if relative.is_absolute() or ".." in relative.parts:
            results.append(("missing", name, "-", "-", "not a path inside the root"))
            continue
        if has_symlink(root, relative):
            results.append(("unknown", name, "-", "-", "symlink"))
            continue
        if not path.is_file():
            results.append(("missing", name, "-", "-", "not a regular file inside the root"))
            continue
        digest = sha256(path)
        rule = exact.get(name)
        if rule is None and not name.startswith(exact_only):
            rule = next((rule for rule in patterns if fnmatch.fnmatchcase(name, rule["pattern"])), None)
        if rule is None:
            results.append(("unknown", name, "-", digest, "no allowlist rule"))
        elif rule.get("sha256", digest) != digest:
            results.append(("mismatch", name, rule["licence"], digest, f"allowlist sha256 {rule['sha256']}"))
        else:
            status = "cleared" if rule["distributable"] else "excluded"
            results.append((status, name, rule["licence"], digest, rule["reason"]))
            if rule["distributable"]:
                for credit in rule.get("credit_files", ()):
                    credits.setdefault(credit, []).append(name)
    cleared = {result[1] for result in results if result[0] == "cleared"}
    for credit, needed_by in credits.items():
        if credit not in cleared:
            results.append(("attribution", credit, "-", "-", f"attribution missing for {', '.join(needed_by)}"))
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=PROJECT, help="repository folder")
    parser.add_argument("--allowlist", type=Path, help="default: ROOT/data/redistribution.json")
    parser.add_argument("--stdin", action="store_true", help="read file paths, one per line, from stdin")
    parser.add_argument("-z", "--null", action="store_true", help="read NUL separated paths from stdin")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    rules, exact_only = load_allowlist(args.allowlist or root / "data/redistribution.json")
    if args.stdin or args.null:
        names = read_names(sys.stdin.buffer, args.null)
    else:
        names = default_files(root)
    results = check(root, names, rules, exact_only)
    for status, name, licence, digest, reason in results:
        print("\t".join((status, licence, digest, name, reason)))
    counts = {status: sum(result[0] == status for result in results)
              for status in ("cleared", "excluded", *FAILURES)}
    print(", ".join(f"{status} {count}" for status, count in counts.items()))
    return 1 if any(counts[status] for status in FAILURES) else 0


if __name__ == "__main__":
    sys.exit(main())
