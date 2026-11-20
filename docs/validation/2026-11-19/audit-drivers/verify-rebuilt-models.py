#!/usr/bin/env python3
"""Audit complete training counts and real classification for selected rebuilt models."""
import argparse
from collections import Counter
import csv
from dataclasses import replace
import gc
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import zipfile

import pandas as pd
from qiime2 import Artifact
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "scripts"))
import build_reference_classifiers as builder
import reference_resources as refs
from taxonomy_provenance import check_taxonomy_provenance


def data_lines(path, suffix):
    with zipfile.ZipFile(path) as archive:
        member, = [name for name in archive.namelist()
                   if name.count("/") == 2 and name.split("/")[1] == "data"
                   and name.endswith(suffix)]
        with archive.open(member) as payload, io.TextIOWrapper(payload) as text:
            yield from text


def sequence_ids(path):
    for line in data_lines(path, ".fasta"):
        if line.startswith(">"):
            yield line[1:].split()[0]


def training_counts(reads, taxonomy):
    rows = csv.reader(data_lines(taxonomy, ".tsv"), delimiter="\t")
    column = next(rows).index("Taxon")
    labels = {}
    for row in rows:
        assert row[0] not in labels, "Duplicate reference taxonomy ID"
        labels[row[0]] = row[column].strip()
    counts, seen = Counter(), set()
    for identifier in sequence_ids(reads):
        assert identifier not in seen, "Duplicate reference sequence ID"
        seen.add(identifier)
        counts[labels[identifier]] += 1
    assert counts, "Empty training references"
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("models", "build-work", "reads", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--classifier", choices=["all", *builder.MODEL_KEYS], default="all")
    args = parser.parse_args()
    selected = list(builder.MODEL_KEYS) if args.classifier == "all" else [args.classifier]
    args.output.mkdir(parents=True)
    query_ids = list(sequence_ids(args.reads))
    assert query_ids and len(query_ids) == len(set(query_ids))
    report = {"status": "incomplete", "environment": builder.environment_versions(),
              "source_revision": subprocess.check_output(
                  ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
              "reads": str(args.reads.resolve()), "requested_classifiers": selected, "models": []}
    report_path = args.output / "model-validation.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    for name in selected:
        started = time.monotonic()
        configured = refs.RESOURCES[builder.MODEL_KEYS[name]]
        model = args.models / configured.filename
        manifest = json.loads(Path(str(model) + ".build.json").read_text())
        identity = manifest["classifier"]
        assert manifest["classifier_id"] == name and identity["filename"] == model.name
        assert identity["sklearn_version"] == refs.CLASSIFIER_SKLEARN_VERSION
        assert builder.recipe_digest(manifest["recipe"], manifest["environment"]) == manifest["recipe_sha256"]
        resource = replace(configured, expected_uuid=identity["uuid"],
                           md5=identity["md5"], sha256=identity["sha256"])
        refs.validate_resource(model, resource)
        workspace = args.build_work / f"{name}-{manifest['recipe_sha256'][:16]}"
        if name == "gg2-2024.09":
            reference_reads, reference_taxonomy = [
                workspace / "sources" / manifest["sources"][item]["filename"]
                for item in ("reference_reads", "reference_taxonomy")]
        else:
            reference_reads, reference_taxonomy = [
                workspace / (item + ".qza") for item in ("reference_reads", "reference_taxonomy")]
        expected = training_counts(reference_reads, reference_taxonomy)
        pipeline = Artifact.load(model).view(Pipeline)
        fitted = pipeline.named_steps["classify"]
        assert len(fitted.classes_) == len(expected) and set(fitted.classes_) == set(expected)
        assert all(float(count) == expected[str(label)]
                   for label, count in zip(fitted.classes_, fitted.class_count_)), "Training class counts differ"
        total = float(fitted.class_count_.sum())
        assert total == sum(expected.values()), "Training sequence totals differ"
        del fitted, pipeline
        gc.collect()  # Release the full model before the classification subprocess loads it.
        taxonomy_path = args.output / f"{name}-taxonomy.qza"
        command = [str(Path(sys.executable).with_name("qiime")), "feature-classifier", "classify-sklearn",
                   "--i-classifier", str(model), "--i-reads", str(args.reads), "--p-n-jobs", "1",
                   "--p-read-orientation", "same", "--p-confidence", "0.7",
                   "--o-classification", str(taxonomy_path)]
        with (args.output / f"{name}-classify.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        artifact = Artifact.load(taxonomy_path)
        artifact.validate(level="min")
        assert str(artifact.type) == "FeatureData[Taxonomy]"
        taxonomy = artifact.view(pd.DataFrame)
        assert {"Taxon", "Confidence"} <= set(taxonomy.columns)
        assert not taxonomy.index.has_duplicates and set(taxonomy.index) == set(query_ids)
        assert all(isinstance(value, str) and value.strip() for value in taxonomy["Taxon"])
        assert all(math.isfinite(float(value)) and 0 <= float(value) <= 1 for value in taxonomy["Confidence"])
        provenance = check_taxonomy_provenance(taxonomy_path, model, reads_path=args.reads, confidence=0.7)
        taxonomy.to_csv(args.output / f"{name}-taxonomy.tsv", sep="\t")
        report["models"].append({"classifier_id": name, "status": "passed", **provenance,
                                 "sha256": identity["sha256"], "training_sequences": int(total),
                                 "training_classes": len(expected), "all_class_counts_match": True,
                                 "reference_reads": str(reference_reads.resolve()),
                                 "reference_taxonomy": str(reference_taxonomy.resolve()),
                                 "features": len(taxonomy), "command": command,
                                 "seconds": round(time.monotonic() - started, 2)})
        report["status"] = "passed" if len(report["models"]) == len(selected) else "incomplete"
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Model validation passed: {report_path}")


if __name__ == "__main__":
    main()
