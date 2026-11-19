#!/usr/bin/env python3
"""Maintainer-only rebuilds of Meta2Data's two fixed legacy classifiers.

Run in the QIIME 2 2026.7 / scikit-learn 1.7.1 environment. This helper is not
called by Meta2Data and accepts no arbitrary reference inputs or classifiers.
SILVA 144 uses its official pretrained classifier and is not rebuilt here.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import zipfile

import reference_resources as refs


SUPPORTED_QIIME = "2026.7"
SCHEMA_VERSION = 1
SOURCE_SILVA_BASE = (
    "https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/")
MODEL_KEYS = {"gg2-2024.09": "greengenes-classifier",
              "silva-138.1": "silva-classifier"}
SILVA_RANKS = ["domain", "phylum", "class", "order", "family", "genus"]

# Preserve the 1.4.2 models' recorded fit parameters except chunk_size:
# smaller training batches reduce peak memory while retaining all references.
# The historical silva-138-99 model uses 138.1 without species labels;
# do not substitute the older 138.0 reference exports.
FIT_PARAMETERS = {
    "classify__alpha": 0.001,
    "classify__chunk_size": 1000,
    "classify__class_prior": "null",
    "classify__fit_prior": False,
    "feat_ext__alternate_sign": False,
    "feat_ext__analyzer": "char_wb",
    "feat_ext__binary": False,
    "feat_ext__decode_error": "strict",
    "feat_ext__encoding": "utf-8",
    "feat_ext__input": "content",
    "feat_ext__lowercase": True,
    "feat_ext__n_features": 8192,
    "feat_ext__ngram_range": "[7, 7]",
    "feat_ext__norm": "l2",
    "feat_ext__preprocessor": "null",
    "feat_ext__stop_words": "null",
    "feat_ext__strip_accents": "null",
    "feat_ext__token_pattern": r"(?u)\b\w\w+\b",
    "feat_ext__tokenizer": "null",
    "verbose": False,
}


def output(name, artifact_id, semantic_type, *, allow_empty_data=False):
    result = {"name": name, "id": artifact_id, "type": semantic_type}
    if allow_empty_data:
        result["allow_empty_data"] = True
    return result


def action(plugin, name, inputs, parameters, outputs):
    return {"plugin": plugin, "action": name, "inputs": inputs,
            "parameters": parameters, "outputs": outputs}


def recipe(classifier_id):
    """Return the fixed, reviewable source and action manifest without I/O."""
    if classifier_id not in MODEL_KEYS:
        raise ValueError(f"No maintainer rebuild recipe for {classifier_id}")
    steps = []
    if classifier_id == "gg2-2024.09":
        sources = {
            "reference_reads": {
                "url": refs.GG2_BASE + "2024.09.backbone.full-length.fna.qza",
                "filename": "2024.09.backbone.full-length.fna.qza",
                "type": "FeatureData[Sequence]",
                "uuid": "5b42d9b6-2f24-4f01-b989-9b4dafca7d5e",
                "payload_md5": "5c5d84ddf524b815ec47a684ace39dcc"},
            "reference_taxonomy": {
                "url": refs.GG2_BASE + "2024.09.backbone.tax.qza",
                "filename": "2024.09.backbone.tax.qza",
                "type": "FeatureData[Taxonomy]",
                "uuid": "b7c3e691-ea51-4547-94dd-f79f49e41a36",
                "payload_md5": "058c1a301903be986cbfb5fdd8bad053"},
        }
        original_uuid = "49ccfb0a-155d-404b-80b4-818d2aeb53b2"
        parse_ranks = None
    else:
        sources = {}
        records = [
            ("rna", "SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz",
             "0b47532bddb5cb2765f38d7420345cbd", "b1a79c7b5a65b7f3b1a4543a8fb3c4e3",
             "FeatureData[RNASequence]", "RNAFASTAFormat"),
            ("taxonomy_tree", "taxonomy/tax_slv_ssu_138.1.tre.gz",
             "287ce923d6f2b18bb4a0e02767c7cb30", "2a232dabdf460da4db143444ae02fdb7",
             "Phylogeny[Rooted]", "NewickFormat"),
            ("taxonomy_ranks", "taxonomy/tax_slv_ssu_138.1.txt.gz",
             "bc864c2147e1661395336a3576ac5062", "65c4258fa6b2fd18feeb7c23a615d104",
             "FeatureData[SILVATaxonomy]", "SILVATaxonomyFormat"),
            ("taxonomy_map", "taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz",
             "245d63d9ef80caf46f50f89d4374200d", "10e46e105a3625d254b1abe90baa8706",
             "FeatureData[SILVATaxidMap]", "SILVATaxidMapFormat"),
        ]
        for source_id, relative, compressed_md5, raw_md5, semantic_type, format_name in records:
            sources[source_id] = {"url": SOURCE_SILVA_BASE + relative,
                                  "filename": Path(relative).name,
                                  "md5": compressed_md5, "uncompressed_md5": raw_md5}
            steps.append({"action": "import", "source": source_id, "format": format_name,
                          "outputs": [output("artifact", source_id, semantic_type)]})
        steps.extend([
            action("rescript", "parse_silva_taxonomy",
                   {name: name for name in ("taxonomy_tree", "taxonomy_map", "taxonomy_ranks")},
                   {"rank_propagation": True, "ranks": SILVA_RANKS,
                    "include_species_labels": False},
                   [output("taxonomy", "taxonomy", "FeatureData[Taxonomy]")]),
            action("rescript", "reverse_transcribe", {"rna_sequences": "rna"}, {},
                   [output("dna_sequences", "dna", "FeatureData[Sequence]")]),
            action("rescript", "cull_seqs", {"sequences": "dna"},
                   {"num_degenerates": 5, "homopolymer_length": 8, "n_jobs": 1},
                   [output("clean_sequences", "clean", "FeatureData[Sequence]")]),
            action("rescript", "filter_seqs_length_by_taxon",
                   {"sequences": "clean", "taxonomy": "taxonomy"},
                   {"labels": ["Archaea", "Bacteria", "Eukaryota"],
                    "min_lens": [900, 1200, 1400], "max_lens": None,
                    "global_min": None, "global_max": None},
                   [output("filtered_seqs", "filtered", "FeatureData[Sequence]"),
                    output("discarded_seqs", "discarded", "FeatureData[Sequence]",
                           allow_empty_data=True)]),
            action("rescript", "dereplicate", {"sequences": "filtered", "taxa": "taxonomy"},
                   {"mode": "uniq", "perc_identity": 1.0, "threads": 1,
                    "rank_handles": SILVA_RANKS + ["species"], "derep_prefix": False},
                   [output("dereplicated_sequences", "reference_reads", "FeatureData[Sequence]"),
                    output("dereplicated_taxa", "reference_taxonomy", "FeatureData[Taxonomy]")]),
        ])
        original_uuid = "70b4b5f4-8fce-40bd-b508-afacbc12a5ed"
        parse_ranks = SILVA_RANKS
    steps.append(action("feature_classifier", "fit_classifier_naive_bayes",
                        {"reference_reads": "reference_reads",
                         "reference_taxonomy": "reference_taxonomy"},
                        dict(FIT_PARAMETERS),
                        [output("classifier", "classifier", "TaxonomicClassifier")]))
    return {"schema_version": SCHEMA_VERSION, "classifier_id": classifier_id,
            "original_classifier_uuid": original_uuid, "taxonomy_parse_ranks": parse_ranks,
            "sklearn_version": refs.CLASSIFIER_SKLEARN_VERSION,
            "sources": sources, "steps": steps}


def environment_versions():
    packages = ("qiime2", "q2-feature-classifier", "rescript", "scikit-learn", "numpy", "scipy")
    versions = {name: importlib.metadata.version(name) for name in packages}
    versions["python"] = platform.python_version()
    for name in ("qiime2", "q2-feature-classifier", "rescript"):
        if not versions[name].startswith(SUPPORTED_QIIME + "."):
            raise RuntimeError(f"{name} {versions[name]} is not the supported {SUPPORTED_QIIME} build")
    if versions["scikit-learn"] != refs.CLASSIFIER_SKLEARN_VERSION:
        raise RuntimeError(f"Activate the scikit-learn {refs.CLASSIFIER_SKLEARN_VERSION} "
                           "Meta2Data environment before building classifiers")
    return versions


def recipe_digest(plan, versions):
    encoded = json.dumps({"recipe": plan, "environment": versions},
                         sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def single_payload_md5(path):
    """Verify that a published input artifact contains the original imported data."""
    with zipfile.ZipFile(path) as archive:
        data = [item for item in archive.infolist()
                if "/data/" in item.filename and not item.is_dir()]
        if len(data) != 1:
            raise ValueError("expected exactly one source artifact data file")
        digest = hashlib.md5()
        with archive.open(data[0]) as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()


def validate_source(path, source):
    if "type" in source:
        refs.validate_artifact(path, source["type"], expected_uuid=source["uuid"])
        if single_payload_md5(path) != source["payload_md5"]:
            raise ValueError("source artifact payload does not match the original classifier input")
    elif refs.file_digest(path, "md5") != source["md5"]:
        raise ValueError("source download does not match its published MD5")


def acquire_source(source, directory, downloader=None):
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / source["filename"]
    if target.exists():
        validate_source(target, source)
    else:
        with tempfile.NamedTemporaryFile(dir=directory, prefix=".download-", delete=False) as tmp:
            staging = Path(tmp.name)
        try:
            (downloader or refs.download)(source["url"], staging)
            validate_source(staging, source)
            os.replace(staging, target)
        finally:
            staging.unlink(missing_ok=True)
    record = {**source, "download_sha256": refs.file_digest(target)}
    if "uncompressed_md5" not in source:
        return target, record
    unpacked = target.with_suffix("")
    if not unpacked.exists():
        with tempfile.NamedTemporaryFile(dir=directory, prefix=".unpack-", delete=False) as tmp:
            staging = Path(tmp.name)
        try:
            with gzip.open(target, "rb") as src, staging.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            if refs.file_digest(staging, "md5") != source["uncompressed_md5"]:
                raise ValueError("uncompressed source differs from the original import provenance")
            os.replace(staging, unpacked)
        finally:
            staging.unlink(missing_ok=True)
    elif refs.file_digest(unpacked, "md5") != source["uncompressed_md5"]:
        raise ValueError(f"Invalid cached source: {unpacked}; remove it and retry")
    record["uncompressed_sha256"] = refs.file_digest(unpacked)
    return unpacked, record


def run_step(step, artifacts, raw_sources, workspace):
    from qiime2 import Artifact

    destinations = {entry["id"]: workspace / (entry["id"] + ".qza")
                    for entry in step["outputs"]}
    def valid(entry):
        kwargs = {"require_nonempty_data": not entry.get("allow_empty_data", False)}
        if entry["type"] == "TaxonomicClassifier":
            kwargs["sklearn_version"] = refs.CLASSIFIER_SKLEARN_VERSION
        return refs.valid_artifact(destinations[entry["id"]], entry["type"], **kwargs)
    if all(valid(entry) for entry in step["outputs"]):
        print(f"Reuse {step['action']}", file=sys.stderr)
        artifacts.update(destinations)
        return
    # Only intermediates in this recipe/environment-specific workspace are
    # replaced. Published models and user analysis outputs are never touched.
    for destination in destinations.values():
        destination.unlink(missing_ok=True)
    print(f"Run {step.get('plugin', 'tools')}.{step['action']}", file=sys.stderr)
    if step["action"] == "import":
        entry = step["outputs"][0]
        imported = Artifact.import_data(entry["type"], str(raw_sources[step["source"]]),
                                        view_type=step["format"])
        imported.save(str(destinations[entry["id"]]))
    else:
        module = importlib.import_module(f"qiime2.plugins.{step['plugin']}.actions")
        method = getattr(module, step["action"])
        inputs = {name: Artifact.load(str(artifacts[artifact_id]))
                  for name, artifact_id in step["inputs"].items()}
        result = method(**inputs, **step["parameters"])
        for entry in step["outputs"]:
            getattr(result, entry["name"]).save(str(destinations[entry["id"]]))
    for entry in step["outputs"]:
        if not valid(entry):
            raise RuntimeError(f"Invalid output from {step['action']}: {entry['id']}")
    artifacts.update(destinations)


def classifier_uuid(path):
    with zipfile.ZipFile(path) as archive:
        roots = {name.split("/")[0] for name in archive.namelist() if "/" in name}
        if len(roots) != 1:
            raise ValueError("classifier must have exactly one artifact root")
        return roots.pop()


def build(classifier_id, output_dir, work_dir, versions):
    plan = recipe(classifier_id)
    fingerprint = recipe_digest(plan, versions)
    workspace = work_dir / f"{classifier_id}-{fingerprint[:16]}"
    workspace.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    resource = refs.RESOURCES[MODEL_KEYS[classifier_id]]
    target = output_dir / resource.filename
    manifest_path = output_dir / (resource.filename + ".build.json")
    if target.exists() or manifest_path.exists():
        raise FileExistsError(f"Output already exists: {target}. Use a new output directory.")
    (workspace / "recipe.json").write_text(json.dumps(plan, indent=2) + "\n")
    artifacts, raw_sources, sources = {}, {}, {}
    for source_id, source in plan["sources"].items():
        path, record = acquire_source(source, workspace / "sources")
        sources[source_id] = record
        (artifacts if "type" in source else raw_sources)[source_id] = path
    for step in plan["steps"]:
        run_step(step, artifacts, raw_sources, workspace)
    model = artifacts["classifier"]
    # A rebuild gets its own UUID and archive hash. Pin those only when the
    # maintainer publishes this build, not to a previously published model.
    manifest = {"schema_version": SCHEMA_VERSION, "classifier_id": classifier_id,
                "built_utc": datetime.now(timezone.utc).isoformat(),
                "recipe_sha256": fingerprint, "environment": versions,
                "classifier": {"filename": resource.filename, "uuid": classifier_uuid(model),
                               "sklearn_version": refs.CLASSIFIER_SKLEARN_VERSION,
                               "sha256": refs.file_digest(model), "md5": refs.file_digest(model, "md5")},
                "sources": sources, "recipe": plan}
    with tempfile.NamedTemporaryFile(dir=output_dir, prefix=".classifier-", delete=False) as tmp:
        staging = Path(tmp.name)
    with tempfile.NamedTemporaryFile(dir=output_dir, prefix=".manifest-", delete=False) as tmp:
        manifest_staging = Path(tmp.name)
    try:
        shutil.copyfile(model, staging)
        if refs.file_digest(staging) != manifest["classifier"]["sha256"]:
            raise RuntimeError("Model changed while copying the finished build")
        manifest_staging.write_text(json.dumps(manifest, indent=2) + "\n")
        os.replace(manifest_staging, manifest_path)
        os.replace(staging, target)
    finally:
        staging.unlink(missing_ok=True)
        manifest_staging.unlink(missing_ok=True)
    print(f"Built {target}")
    print(f"Manifest {manifest_path}")
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classifier", choices=[*MODEL_KEYS, "all"], required=True)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Fixed built-in resource destination, e.g. Meta2Data/references")
    parser.add_argument("--work-dir", type=Path,
                        help="Build intermediates; defaults to <output-dir>/build-work")
    parser.add_argument("--plan", action="store_true",
                        help="Print the fixed recipes without importing QIIME, downloading, or training")
    args = parser.parse_args()
    selected = list(MODEL_KEYS) if args.classifier == "all" else [args.classifier]
    if args.plan:
        print(json.dumps([recipe(name) for name in selected], indent=2))
        return
    try:
        versions = environment_versions()
        destination = args.output_dir.resolve()
        workspace = (args.work_dir or destination / "build-work").resolve()
        for name in selected:
            build(name, destination, workspace, versions)
    except (OSError, ValueError, RuntimeError, importlib.metadata.PackageNotFoundError) as error:
        parser.exit(1, f"Build failed: {error}\n")


if __name__ == "__main__":
    main()
