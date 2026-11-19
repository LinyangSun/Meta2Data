#!/usr/bin/env python3
"""Check taxonomy cache provenance without extracting or loading model payloads."""
import argparse
import math
from pathlib import Path
import sys
import uuid
import zipfile
import zlib

import yaml


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_yaml(archive, name):
    # Metadata/action records are small. Reject a malformed archive claiming
    # a model-sized record rather than inflating it into memory.
    require(archive.getinfo(name).file_size <= 1024 * 1024,
            f"Provenance record is unexpectedly large: {name}")
    # QIIME's !ref and !cite scalars stay text, never Python objects.
    return yaml.load(archive.read(name), Loader=yaml.BaseLoader)


def valid_uuid(value):
    try:
        return isinstance(value, str) and str(uuid.UUID(value)) == value
    except ValueError:
        return False


def artifact_uuid(archive, path, semantic_type):
    names = [name for name in archive.namelist()
             if name.count("/") == 1 and name.endswith("/metadata.yaml")]
    require(len(names) == 1, f"Expected one artifact metadata record: {path}")
    metadata = read_yaml(archive, names[0])
    require(isinstance(metadata, dict), f"Invalid artifact metadata: {path}")
    value = metadata.get("uuid")
    require(valid_uuid(value) and names[0].split("/")[0] == value,
            f"Artifact UUID and archive root differ or are invalid: {path}")
    require(metadata.get("type") == semantic_type,
            f"{path}: expected {semantic_type}, got {metadata.get('type')}")
    return value


def one_value(action, section, name):
    records = action.get(section)
    require(isinstance(records, list), f"Taxonomy provenance has no {section} records")
    values = [entry[name] for entry in records if isinstance(entry, dict) and name in entry]
    require(len(values) == 1, f"Taxonomy provenance must contain one {name} record, found {values}")
    return values[0]


def check_taxonomy_provenance(taxonomy_path, classifier_path, *, reads_path=None, confidence=None):
    """Match the producing action to the requested model, reads and threshold.

    Reference preparation separately verifies published model identity. Optional
    reads/confidence arguments allow inspection of already validated results;
    production cache reuse always supplies both.
    """
    try:
        with zipfile.ZipFile(classifier_path) as archive:
            classifier = artifact_uuid(archive, classifier_path, "TaxonomicClassifier")
        reads = None
        if reads_path is not None:
            with zipfile.ZipFile(reads_path) as archive:
                reads = artifact_uuid(archive, reads_path, "FeatureData[Sequence]")
        with zipfile.ZipFile(taxonomy_path) as archive:
            taxonomy = artifact_uuid(archive, taxonomy_path, "FeatureData[Taxonomy]")
            provenance = read_yaml(archive, f"{taxonomy}/provenance/action/action.yaml")
    except (OSError, KeyError, EOFError, zipfile.BadZipFile, zlib.error, yaml.YAMLError,
            UnicodeError, RuntimeError, NotImplementedError) as error:
        raise ValueError(f"Cannot read classifier provenance for {taxonomy_path}: {error}") from error
    action = provenance.get("action") if isinstance(provenance, dict) else None
    require(isinstance(action, dict) and action.get("type") == "method" and
            action.get("plugin") == "environment:plugins:feature-classifier" and
            action.get("action") == "classify_sklearn",
            f"Taxonomy was not produced directly by feature-classifier.classify_sklearn: {taxonomy_path}")
    observed_classifier = one_value(action, "inputs", "classifier")
    require(observed_classifier == classifier,
            f"Taxonomy classifier UUID mismatch: expected {classifier} from {classifier_path}, "
            f"found {observed_classifier} in {taxonomy_path}")
    observed_reads = one_value(action, "inputs", "reads")
    require(valid_uuid(observed_reads), f"Invalid taxonomy reads UUID: {observed_reads!r}")
    if reads is not None:
        require(observed_reads == reads,
                f"Taxonomy reads UUID mismatch: expected {reads}, found {observed_reads}")
    observed_confidence = one_value(action, "parameters", "confidence")
    try:
        observed_confidence = float(observed_confidence)
        expected_confidence = None if confidence is None else float(confidence)
    except (TypeError, ValueError) as error:
        raise ValueError("Taxonomy confidence must be a finite numeric value between 0 and 1") from error
    require(math.isfinite(observed_confidence) and 0 <= observed_confidence <= 1,
            f"Invalid taxonomy confidence: {observed_confidence}")
    if expected_confidence is not None:
        require(math.isfinite(expected_confidence) and 0 <= expected_confidence <= 1,
                f"Invalid requested confidence: {expected_confidence}")
        require(observed_confidence == expected_confidence,
                f"Taxonomy confidence mismatch: expected {expected_confidence}, found {observed_confidence}")
    return {"taxonomy_uuid": taxonomy, "classifier_uuid": classifier,
            "reads_uuid": observed_reads, "confidence": observed_confidence}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--classifier", required=True, type=Path)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--confidence", required=True, type=float)
    args = parser.parse_args(argv)
    try:
        check_taxonomy_provenance(args.taxonomy, args.classifier,
                                  reads_path=args.reads, confidence=args.confidence)
    except ValueError as error:
        print(f"Cached taxonomy rejected: {error}; recomputing.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
