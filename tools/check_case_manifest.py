"""Check case manifest provenance and split independence using committed metadata."""

import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "data/case_manifest.json"
INVENTORY = PROJECT / "data/manifest.csv"
SPLITS = {"development", "final_holdout", "excluded"}


def _index(rows, kind):
    result = {}
    for row in rows:
        identity = row["id"]
        if not isinstance(identity, str) or not identity:
            raise ValueError(f"Invalid {kind} id: {identity!r}")
        if identity in result:
            raise ValueError(f"Duplicate {kind} id: {identity}")
        result[identity] = row
    return result


def _file_key(inventory, key, field, digest=None):
    key = tuple(key)
    if key not in inventory:
        raise ValueError(f"{field} missing from inventory: {key}")
    if not inventory[key] or (digest is not None and inventory[key] != digest):
        raise ValueError(f"{field} sha256 mismatch: {key}")
    return key


def _read_inventory(path):
    contents = path.read_bytes()
    if path == INVENTORY:
        rules = json.loads((PROJECT / "data/redistribution.json").read_text(encoding="utf-8"))["rules"]
        entry = next((rule for rule in rules if rule.get("path") == "data/manifest.csv"), {})
        if hashlib.sha256(contents).hexdigest() != entry.get("sha256"):
            raise ValueError("Committed inventory sha256 mismatch")
    inventory = {}
    for row in csv.DictReader(io.StringIO(contents.decode("utf-8"), newline="")):
        key = row["deposit_id"], row["archive"], row["path"]
        if key in inventory:
            raise ValueError(f"Duplicate inventory key: {key}")
        inventory[key] = row["sha256"]
    return inventory


def validate(manifest_path, inventory_csv, *, reference_inventory=INVENTORY):
    """Return sorted components, checking processed digests against the reference inventory."""
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    inventory_path = Path(inventory_csv).resolve()
    reference_path = Path(reference_inventory).resolve()
    inventory = _read_inventory(inventory_path)
    reference = inventory if inventory_path == reference_path else _read_inventory(reference_path)

    runs = _index(manifest["runs"], "run")
    cases = _index(manifest["cases"], "case")
    preparations = _index(manifest["preparations"], "preparation")
    edges = {(kind, identity): set() for kind, rows in
             (("run", runs), ("case", cases), ("preparation", preparations)) for identity in rows}

    def connect(left, right):
        if right not in edges:
            raise ValueError(f"Dangling {right[0]} reference {right[1]!r} from {left[1]}")
        edges[left].add(right)
        edges[right].add(left)

    files, digests, acquisitions, fixtures = {}, {}, {}, {}
    for identity, run in runs.items():
        node = "run", identity
        digest = run["sha256"]
        if not isinstance(digest, str) or not digest:
            raise ValueError(f"Invalid raw_key sha256 for {identity}")
        raw = _file_key(inventory, run["raw_key"], "raw_key", digest)
        connect(node, files.setdefault(raw, node))
        connect(node, digests.setdefault(digest, node))
        if run["processed_key"] is not None:
            sibling = _file_key(reference, run["processed_key"], "processed_key", run.get("processed_sha256"))
            _file_key(inventory, sibling, "processed_key", reference[sibling])
            connect(node, files.setdefault(sibling, node))
        acquisition = run["acquisition_id"]
        if acquisition is not None:
            if isinstance(acquisition, dict):
                acquisition = tuple(sorted(acquisition.items()))
            connect(node, acquisitions.setdefault((raw[0], acquisition), node))
        if run["alias_of"] is not None:
            connect(node, ("run", run["alias_of"]))
        if run["preparation_id"] is not None:
            connect(node, ("preparation", run["preparation_id"]))
        for path in run["fixture_aliases"]:
            owners = fixtures.setdefault(path, set())
            if owners:
                connect(node, next(iter(owners)))
            owners.add(node)

    for link in manifest["links"]:
        left, right = (("preparation", identity) for identity in link["preparation_ids"])
        if left not in edges:
            raise ValueError(f"Dangling preparation reference {left[1]!r} in link")
        connect(left, right)

    for identity, case in cases.items():
        if case["split"] not in SPLITS:
            raise ValueError(f"Unknown split {case['split']!r} for {identity}")
        if not isinstance(case["run_ids"], list) or not case["run_ids"]:
            raise ValueError(f"Case {identity} has no run references")
        if case.get("kind") == "window":
            path = case["path"]
            if not isinstance(path, str):
                raise ValueError(f"Invalid window path for {identity}")
            if path.startswith("tests/fixtures/") and not any(
                    ("run", run_id) in fixtures.get(path, ()) for run_id in case["run_ids"]):
                raise ValueError(f"Fixture window {identity} has no owner run for {path}")
        for run_id in case["run_ids"]:
            connect(("case", identity), ("run", run_id))

    for identity, pair in _index(manifest["pilot_pairs"], "pilot pair").items():
        ordered = []
        for role in ("sample", "ec"):
            member = pair.get(role)
            run = runs.get(member.get("run_id")) if isinstance(member, dict) else None
            if run is None or run["raw_key"] != member.get("raw_key"):
                raise ValueError(f"Invalid {role} run reference in pilot pair {identity}")
            ordered.append(run["id"])
        if identity not in cases or cases[identity]["run_ids"] != ordered:
            raise ValueError(f"Pilot pair {identity} must match its case's ordered run_ids")

    exposed, unresolved = set(), []
    for exposure in manifest["exposures"]:
        if not isinstance(exposure, dict):
            raise ValueError("Exposure must be an object")
        run_ids = exposure.get("run_ids")
        if not run_ids:
            unresolved.append(exposure["scope"])
            continue
        if not isinstance(run_ids, list) or any(run_id not in runs for run_id in run_ids):
            raise ValueError(f"Dangling exposure run reference: {exposure['scope']}")
        exposed.update(run_ids)

    components, seen = [], set()
    for node in sorted(edges):
        if node in seen:
            continue
        pending, component = [node], []
        seen.add(node)
        while pending:
            member = pending.pop()
            component.append(member)
            for neighbor in edges[member] - seen:
                seen.add(neighbor)
                pending.append(neighbor)
        component = tuple(sorted(component))
        components.append(component)
        splits = {cases[identity]["split"] for kind, identity in component if kind == "case"}
        if "final_holdout" not in splits:
            continue
        label = next(identity for kind, identity in component
                     if kind == "case" and cases[identity]["split"] == "final_holdout")
        if "development" in splits:
            raise ValueError(f"Split conflict for {label}: component contains development and final_holdout")
        run_ids = [identity for kind, identity in component if kind == "run"]
        if any(runs[identity]["preparation_id"] is None for identity in run_ids):
            raise ValueError(f"Final holdout {label} has a null preparation_id in its component")
        if unresolved:
            raise ValueError(f"Final holdout {label} blocked by unresolved exposure: {unresolved[0]}")
        if exposed.intersection(run_ids):
            raise ValueError(f"Final holdout {label} has an exposed run in its component")
    return tuple(components)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--inventory", type=Path, default=INVENTORY)
    args = parser.parse_args()
    try:
        components = validate(args.manifest, args.inventory)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        detail = " ".join(str(error).splitlines())
        print(f"case manifest: INVALID ({detail})", file=sys.stderr)
        return 1
    print(f"case manifest: OK ({len(components)} components, {sum(map(len, components))} nodes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
