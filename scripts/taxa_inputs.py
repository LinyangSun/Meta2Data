#!/usr/bin/env python3
"""Collect paired QIIME artifacts and invalidate dependent TAXA caches."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import uuid
import zipfile

RESULT = re.compile(r"^(.+)-(dada2|vsearch)-final-(table|rep-seqs)\.qza$")
EXPECTED_TYPES = {"table": "FeatureTable[Frequency]", "rep_seqs": "FeatureData[Sequence]"}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def artifact(path, kind):
    """Check archive integrity and semantic type without importing QIIME."""
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        metadata_paths = [name for name in archive.namelist()
                          if len(name.split("/")) == 2 and name.endswith("/metadata.yaml")]
        if len(metadata_paths) != 1:
            raise ValueError("expected one root metadata.yaml")
        metadata = {}
        for line in archive.read(metadata_paths[0]).decode().splitlines():
            if ":" in line and not line.startswith(" "):
                key, value = line.split(":", 1)
                metadata[key] = value.strip().strip("'\"")
        artifact_uuid = str(uuid.UUID(metadata.get("uuid", "")))
        if artifact_uuid != metadata_paths[0].split("/")[0]:
            raise ValueError("artifact UUID does not match archive root")
        if metadata.get("type") != EXPECTED_TYPES[kind]:
            raise ValueError(f"expected {EXPECTED_TYPES[kind]}, found {metadata.get('type')}")
        bad = archive.testzip()
        if bad:
            raise ValueError(f"archive integrity check failed: {bad}")
        if not any(name.startswith(f"{artifact_uuid}/data/") and not name.endswith("/")
                   for name in archive.namelist()):
            raise ValueError("artifact contains no data files")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.absolute()), "realpath": str(path.resolve()),
            "uuid": artifact_uuid, "sha256": digest.hexdigest()}


def discover(directory):
    pairs = {}
    visited = set()
    for root, directories, files in os.walk(directory, followlinks=True):
        realroot = str(Path(root).resolve())
        if realroot in visited or (Path(root) / "taxa-run-state.json").is_file():
            directories[:] = []
            continue
        visited.add(realroot)
        # A final-* prefix is also a valid local dataset name. Only a TAXA
        # state marker identifies an owned aggregate output directory.
        directories[:] = sorted(name for name in directories
                                if name not in {"tmp", ".tmpdir"})
        for name in sorted(files):
            match = RESULT.fullmatch(name)
            if match:
                dataset, method, kind = match.groups()
                key = (str(Path(root)), dataset, method)
                pairs.setdefault(key, {})[kind.replace("-", "_")] = Path(root) / name
    return pairs


def collect(directory, method=None):
    directory = Path(directory)
    audit = {"input_directory": str(directory.resolve()), "method": method,
             "datasets": [], "duplicates": [], "skipped": [], "errors": []}
    if not directory.is_dir():
        audit["errors"].append(f"Input directory does not exist: {directory}")
        return audit
    pairs = discover(directory)
    complete = {}
    for key, paths in pairs.items():
        if len(paths) != 2:
            audit["skipped"].append({"directory": key[0], "dataset": key[1], "method": key[2],
                                     "reason": "Missing matching table or representative-sequences artifact"})
        else:
            complete[key] = paths
    methods = {key[2] for key in complete}
    if not method and len(methods) > 1:
        audit["errors"].append("Both dada2 and vsearch results are present. Select --dada2 or --vsearch.")
        return audit
    if not method and methods:
        method = next(iter(methods))
    audit["method"] = method
    seen_uuids, seen_pairs, seen_datasets, seen_members = {}, {}, {}, {}
    for (folder, dataset, found_method), paths in sorted(complete.items()):
        if found_method != method:
            audit["skipped"].append({"directory": folder, "dataset": dataset, "method": found_method,
                                     "reason": "Different feature method"})
            continue
        try:
            record = {"dataset": dataset, "method": method,
                      **{kind: artifact(path, kind) for kind, path in paths.items()}}
        except (ValueError, OSError, zipfile.BadZipFile, RuntimeError) as error:
            audit["skipped"].append({"directory": folder, "dataset": dataset, "method": method,
                                     "reason": f"Invalid artifact: {error}"})
            continue
        pair_key = tuple(record[kind]["uuid"] for kind in EXPECTED_TYPES)
        pair_hash = tuple(record[kind]["sha256"] for kind in EXPECTED_TYPES)
        conflict = False
        for kind in EXPECTED_TYPES:
            info = record[kind]
            previous = seen_uuids.get(info["uuid"])
            if previous and previous["sha256"] != info["sha256"]:
                audit["errors"].append(f"Conflicting content for artifact UUID {info['uuid']}: "
                                       f"{previous['path']} and {info['path']}")
                conflict = True
            seen_uuids[info["uuid"]] = info
            paired_with = seen_members.get(info["uuid"])
            if paired_with and paired_with != pair_key:
                audit["errors"].append(f"Conflicting pairing for artifact UUID {info['uuid']}: "
                                       f"{info['path']} reuses only one member of an existing result pair.")
                conflict = True
        if dataset in seen_datasets and seen_datasets[dataset] != pair_hash:
            audit["errors"].append(f"Dataset name '{dataset}' has different results in multiple folders; "
                                   "use distinct dataset names or select a narrower input directory.")
            conflict = True
        seen_datasets[dataset] = pair_hash
        if conflict:
            continue
        if pair_key in seen_pairs:
            audit["duplicates"].append({"dataset": dataset, "directory": folder,
                                         "same_as": seen_pairs[pair_key], "reason": "Same artifact UUID pair"})
            continue
        seen_pairs[pair_key] = str(paths["table"].parent)
        for artifact_uuid in pair_key:
            seen_members[artifact_uuid] = pair_key
        audit["datasets"].append(record)
    if not audit["datasets"] and not audit["errors"]:
        audit["errors"].append("No valid complete result pairs found. Expected "
                               "<id>-dada2-final-table.qza and <id>-dada2-final-rep-seqs.qza "
                               "(or the corresponding vsearch pair).")
    return audit


def database_signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    # Reference databases can be many GB; do not hash them on every restart.
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


CACHE_OUTPUTS = {
    "merge": ["tmp/mergedTable.qza", "tmp/mergedRepSeqs.qza", "tmp/mergedTableSummary.qzv"],
    "orient": ["orientedTable.qza", "orientedRepSeqs.qza", "tmp/orientedTable.qza",
               "tmp/orientedRepSeqs.qza", "tmp/orientedTableSummary.qzv", "tmp/unmatchedRepSeqs.qza"],
    "taxonomy": ["gg2Taxonomy.qza", "silvaTaxonomy.qza"],
    "tree": ["seppTree.qza", "tmp/seppPlacements.qza", "treeFilteredTable.qza",
             "treeFilteredRepSeqs.qza", "tmp/treeUnplacedTable.qza", "tmp/treeFilteredTableSummary.qzv",
             "tmp/denovoAlignment.qza", "tmp/denovoMaskedAlignment.qza", "tmp/denovoUnrootedTree.qza",
             "denovoRootedTree.qza"],
}


def prepare(collection, final_dir, orient_ref, classifier, confidence, single_v=False,
            sepp_ref=None, db_label="gg2", notree=False):
    if single_v and notree:
        raise ValueError("--single-v and --notree are mutually exclusive")
    tree_mode = "none" if notree else "denovo" if single_v else "sepp"
    if db_label not in {"gg2", "silva"}:
        raise ValueError(f"Unknown taxonomy database label: {db_label}")
    audit = json.loads(Path(collection).read_text())
    if audit["errors"] or not audit["datasets"]:
        raise ValueError("Cannot prepare TAXA from an unsuccessful collection")
    final_dir = Path(final_dir)
    final_dir.mkdir(parents=True, exist_ok=True)
    inputs = sorted((row["dataset"], row["table"]["uuid"], row["table"]["sha256"],
                     row["rep_seqs"]["uuid"], row["rep_seqs"]["sha256"])
                    for row in audit["datasets"])
    merge = fingerprint({"method": audit["method"], "inputs": inputs})
    orient = fingerprint({"merge": merge, "reference": database_signature(orient_ref)})
    if tree_mode == "sepp" and not sepp_ref:
        raise ValueError("SEPP runs require a SEPP reference")
    classifier_info = database_signature(classifier)
    dependencies = {
        "merge": merge,
        "orient": orient,
        "taxonomy": fingerprint({"orient": orient, "classifier": classifier_info,
                                  "confidence": confidence}),
        "tree": fingerprint({"orient": orient, "tree_mode": tree_mode,
                              "reference": database_signature(sepp_ref) if tree_mode == "sepp" else None}),
    }
    state_path = final_dir / "taxa-run-state.json"
    previous = json.loads(state_path.read_text()) if state_path.exists() else {}
    previous_dependencies = previous.get("dependencies", {})
    taxonomy_dependencies = previous.get("taxonomy_dependencies", {})
    if previous_dependencies.get("orient") != orient:
        # Both classifiers depend on the same oriented feature sequences.
        for relative in CACHE_OUTPUTS["taxonomy"]:
            (final_dir / relative).unlink(missing_ok=True)
        taxonomy_dependencies = {}
    invalidated = []
    for stage in CACHE_OUTPUTS:
        if stage == "taxonomy":
            cached = taxonomy_dependencies.get(db_label, {}).get("fingerprint")
            outputs = [f"{db_label}Taxonomy.qza"]
        else:
            cached = previous_dependencies.get(stage)
            outputs = CACHE_OUTPUTS[stage]
        if cached != dependencies[stage]:
            invalidated.append(stage)
            for relative in outputs:
                (final_dir / relative).unlink(missing_ok=True)
    taxonomy_dependencies[db_label] = {"fingerprint": dependencies["taxonomy"],
                                       "confidence": confidence, "classifier": classifier_info}
    shutil.copyfile(collection, final_dir / "collection.json")
    write_json(state_path, {"method": audit["method"], "singleV": single_v,
                            "notree": notree, "tree_mode": tree_mode, "confidence": confidence,
                            "dependencies": dependencies,
                            "taxonomy_dependencies": taxonomy_dependencies})
    return invalidated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    collector = commands.add_parser("collect")
    collector.add_argument("--input", required=True)
    collector.add_argument("--method", choices=["dada2", "vsearch"])
    collector.add_argument("--output", required=True)
    preparer = commands.add_parser("prepare")
    preparer.add_argument("--collection", required=True)
    preparer.add_argument("--final-dir", required=True)
    preparer.add_argument("--orient-ref", required=True)
    preparer.add_argument("--classifier", required=True)
    preparer.add_argument("--confidence", type=float, required=True)
    preparer.add_argument("--db-label", choices=["gg2", "silva"], default="gg2")
    tree_options = preparer.add_mutually_exclusive_group()
    tree_options.add_argument("--single-v", action="store_true")
    tree_options.add_argument("--notree", "--no-tree", action="store_true")
    preparer.add_argument("--sepp-ref")
    paths = commands.add_parser("paths")
    paths.add_argument("--collection", required=True)
    paths.add_argument("--kind", choices=EXPECTED_TYPES, required=True)
    args = parser.parse_args()
    try:
        if args.action == "collect":
            audit = collect(args.input, args.method)
            write_json(args.output, audit)
            for row in audit["skipped"]:
                print(f"Skipping {row['dataset']}: {row['reason']}", file=sys.stderr)
            for error in audit["errors"]:
                print(f"Error: {error}", file=sys.stderr)
            print(f"Collection audit: {args.output}", file=sys.stderr)
            if audit["errors"]:
                return 1
            print(f"Collected {len(audit['datasets'])} {audit['method']} datasets; "
                  f"removed {len(audit['duplicates'])} duplicate pairs.", file=sys.stderr)
            print(audit["method"])
        elif args.action == "prepare":
            invalidated = prepare(args.collection, args.final_dir, args.orient_ref, args.classifier,
                                  args.confidence, args.single_v, args.sepp_ref, args.db_label, args.notree)
            if invalidated:
                print("Recompute affected TAXA steps: " + ", ".join(invalidated))
        else:
            audit = json.loads(Path(args.collection).read_text())
            for row in audit["datasets"]:
                sys.stdout.buffer.write(os.fsencode(row[args.kind]["path"]) + b"\0")
    except (ValueError, OSError, zipfile.BadZipFile) as error:
        parser.exit(2, f"TAXA input error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
