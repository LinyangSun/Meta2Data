#!/usr/bin/env python3
"""Validate the actual Meta2Data runtime before an environment or SIF is released.

Default: all unit tests and a small real workflow matrix. --smoke checks a
packaged image with PE DADA2/VSEARCH, every built-in classifier, and both trees.
--checks-only omits fixtures; --environment-only also omits unit tests.
Missing tools, unavailable classifiers, skipped unit tests, and missing outputs
are failures. Synthetic inputs test compatibility, not taxonomic accuracy.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
import uuid
import zipfile

ROOT = Path(__file__).resolve().parents[2]
CLASSIFIERS = {
    "gg2-2024.09": ("greengenes-classifier", "gg2"),
    "silva-138.1": ("silva-classifier", "silva"),
    "silva-144": ("silva144-classifier", "silva144"),
}
VERSIONS = {
    "qiime2": "2026.7.0", "rachis": "2026.7.0", "q2cli": "2026.7.0",
    "q2-dada2": "2026.7.0", "q2-demux": "2026.7.0",
    "q2-feature-classifier": "2026.7.0", "q2-feature-table": "2026.7.0",
    "q2-fragment-insertion": "2026.7.0", "q2-phylogeny": "2026.7.0",
    "q2-quality-filter": "2026.7.0",
    "rescript": "2026.7.0", "scikit-learn": "1.7.1",
}
ACTIONS = {
    "dada2": ("denoise_paired", "denoise_single", "denoise_pyro", "denoise_ccs"),
    "demux": ("summarize",),
    "feature-table": ("merge", "merge_seqs", "filter_features", "filter_seqs", "summarize"),
    "feature-classifier": ("classify_sklearn",),
    "fragment-insertion": ("sepp", "filter_features"),
    "phylogeny": ("align_to_tree_mafft_fasttree",),
    "quality-filter": ("q_score",),
    "rescript": ("orient_seqs",),
}
BIO_BINARIES = ("qiime", "python", "python3", "fastp", "cutadapt", "vsearch",
                "blastn", "makeblastdb", "samtools", "chopper", "minimap2",
                "racon", "prefetch", "fasterq-dump")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_identifiers(table_ids, sequence_ids, taxonomy_ids=None, tree_ids=None):
    """Compare actual IDs, so equal counts cannot conceal different features."""
    table_ids, sequence_ids = set(table_ids), set(sequence_ids)
    require(bool(table_ids), "No features survived the compatibility fixture")
    require(table_ids == sequence_ids, "Feature table and representative-sequence IDs differ")
    if taxonomy_ids is not None:
        require(table_ids <= set(taxonomy_ids), "Features are missing taxonomy records")
    if tree_ids is not None:
        require(table_ids <= set(tree_ids), "Features are missing from the final tree")


def _artifact_metadata(archive, path, semantic_type):
    """Read identity without extracting a potentially large classifier payload."""
    import yaml
    names = [name for name in archive.namelist()
             if name.count("/") == 1 and name.endswith("/metadata.yaml")]
    require(len(names) == 1, f"Expected one artifact metadata record: {path}")
    metadata = yaml.load(archive.read(names[0]), Loader=yaml.BaseLoader)
    require(isinstance(metadata, dict), f"Invalid artifact metadata: {path}")
    artifact_uuid = metadata.get("uuid")
    try:
        valid_uuid = isinstance(artifact_uuid, str) and str(uuid.UUID(artifact_uuid)) == artifact_uuid
    except ValueError:
        valid_uuid = False
    require(valid_uuid and names[0].split("/")[0] == artifact_uuid,
            f"Artifact UUID and archive root differ or are invalid: {path}")
    require(metadata.get("type") == semantic_type,
            f"{path}: expected {semantic_type}, got {metadata.get('type')}")
    return artifact_uuid


def check_taxonomy_classifier(taxonomy_path, classifier_path):
    """Require the output action to reference the selected classifier's UUID.

    Resource preparation independently checks published classifier identity.
    This check catches a pipeline using another model despite the expected
    output filename; it does not replace missing release UUID/hash pins.
    """
    import yaml
    try:
        with zipfile.ZipFile(classifier_path) as archive:
            classifier_uuid = _artifact_metadata(archive, classifier_path, "TaxonomicClassifier")
        with zipfile.ZipFile(taxonomy_path) as archive:
            taxonomy_uuid = _artifact_metadata(archive, taxonomy_path, "FeatureData[Taxonomy]")
            action = yaml.load(archive.read(
                f"{taxonomy_uuid}/provenance/action/action.yaml"), Loader=yaml.BaseLoader)
    except (OSError, KeyError, zipfile.BadZipFile, yaml.YAMLError) as error:
        raise RuntimeError(f"Cannot read classifier provenance for {taxonomy_path}: {error}") from error
    action = action.get("action") if isinstance(action, dict) else None
    require(isinstance(action, dict) and action.get("type") == "method" and
            action.get("plugin") == "environment:plugins:feature-classifier" and
            action.get("action") == "classify_sklearn",
            f"Taxonomy was not produced directly by feature-classifier.classify_sklearn: {taxonomy_path}")
    inputs = action.get("inputs")
    require(isinstance(inputs, list), f"Taxonomy provenance has no input records: {taxonomy_path}")
    classifier_inputs = [entry["classifier"] for entry in inputs
                         if isinstance(entry, dict) and "classifier" in entry]
    require(classifier_inputs == [classifier_uuid],
            f"Taxonomy classifier UUID mismatch: expected {classifier_uuid} from {classifier_path}, "
            f"found {classifier_inputs} in {taxonomy_path}")
    return {"taxonomy_uuid": taxonomy_uuid, "classifier_uuid": classifier_uuid}


def check_silva144_ranks(taxa):
    """Check the observed SILVA 144 schema without requiring full-depth calls."""
    prefixes = ("d__", "k__", "p__", "c__", "o__", "f__", "g__")
    rank_counts, unassigned = set(), 0
    for index, taxon in enumerate(taxa):
        require(isinstance(taxon, str) and bool(taxon.strip()),
                f"SILVA 144 taxonomy row {index} is empty or not text: {taxon!r}")
        taxon = taxon.strip()
        if taxon == "Unassigned":
            rank_counts.add(0)
            unassigned += 1
            continue
        lineage = [rank.strip() for rank in taxon.split(";")]
        require(len(lineage) <= len(prefixes) and
                all(rank.startswith(prefix) for rank, prefix in zip(lineage, prefixes)),
                f"SILVA 144 rank order must be d/k/p/c/o/f/g, possibly truncated: {taxon!r}")
        rank_counts.add(len(lineage))
    return {"observed_rank_counts": sorted(rank_counts), "unassigned_features": unassigned,
            "expected_rank_prefixes": list(prefixes)}


class Validation:
    def __init__(self, work, mode, threads):
        self.work, self.mode, self.threads = work, mode, threads
        require(not work.exists() or (work.is_dir() and not any(work.iterdir())),
                "--work-dir must be absent or empty; previous evidence is never overwritten")
        work.mkdir(parents=True, exist_ok=True)
        (work / "logs").mkdir()
        self.env = dict(os.environ, PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1")
        self.report = {
            "schema_version": 1, "mode": mode, "status": "running",
            "started": datetime.now(timezone.utc).isoformat(),
            "source": str(ROOT), "python_executable": sys.executable,
            "env_yml_sha256": hashlib.sha256((ROOT / "env.yml").read_bytes()).hexdigest(),
            "environment_lock_sha256": hashlib.sha256(
                (ROOT / "environment-linux-64.lock").read_bytes()).hexdigest(),
            "steps": [], "products": {},
        }
        self.save()

    def save(self):
        (self.work / "validation.json").write_text(
            json.dumps(self.report, indent=2, sort_keys=True) + "\n")

    def run(self, name, command, *, cwd=None, env=None):
        command = [str(value) for value in command]
        log = self.work / "logs" / (name + ".log")
        print(f"[{name}] {shlex.join(command)}", flush=True)
        started = time.monotonic()
        step = {"name": name, "command": command, "log": str(log)}
        self.report["steps"].append(step)
        self.save()
        with log.open("w") as stream:
            result = subprocess.run(command, cwd=cwd or ROOT, env=env or self.env,
                                    stdout=stream, stderr=subprocess.STDOUT)
        step.update(returncode=result.returncode, seconds=round(time.monotonic() - started, 2))
        self.save()
        if result.returncode:
            tail = "\n".join(log.read_text(errors="replace").splitlines()[-35:])
            raise RuntimeError(f"{name} failed (exit {result.returncode}); see {log}\n{tail}")
        return log

    def environment(self):
        require(sys.version_info[:2] == (3, 12), "QIIME 2 2026.7 requires this project's Python 3.12 environment")
        found = {name: version(name) for name in VERSIONS}
        self.report["versions"] = found
        require(found == VERSIONS, f"Runtime versions do not match env.yml: {found}")
        for module in ("Bio", "pandas", "numpy", "biom", "openpyxl", "requests", "psutil"):
            importlib.import_module(module)
        prefix = Path(sys.prefix).resolve()
        binaries = {}
        for name in (*BIO_BINARIES, "curl", "wget", "unzip", "gzip", "bash"):
            executable = shutil.which(name)
            require(executable, f"Required executable is missing: {name}")
            binaries[name] = executable
            if name in BIO_BINARIES:
                require(Path(executable).resolve().is_relative_to(prefix),
                        f"{name} resolves outside the active runtime: {executable}")
        installed = shutil.which("Meta2Data")
        require(installed and Path(installed).resolve() == (ROOT / "bin/Meta2Data").resolve(),
                "Meta2Data on PATH must be the entry point from the source or SIF being validated")
        self.report["binaries"] = binaries

        sys.path.insert(0, str(ROOT / "scripts"))
        import install_sra_toolkit as sra
        manifest = prefix / "opt" / sra.DIRECTORY / sra.MANIFEST
        require(manifest.is_file(), f"SRA Toolkit installation record is missing: {manifest}")
        sra_record = json.loads(manifest.read_text())
        require(sra_record.get("version") == sra.VERSION and
                sra_record.get("archive", {}).get("sha256") == sra.ARCHIVE_SHA256,
                "SRA Toolkit installation record does not match the pinned component")
        self.report["sra_toolkit"] = sra_record

        from qiime2.sdk import PluginManager
        plugins = {plugin.name.replace("_", "-"): plugin for plugin in PluginManager().plugins.values()}
        for name, actions in ACTIONS.items():
            require(name in plugins, f"Missing QIIME plugin: {name}")
            for action in actions:
                require(action in plugins[name].actions, f"Missing QIIME action: {name}.{action}")
        for action in ACTIONS["dada2"]:
            require("base_transition_stats" in plugins["dada2"].actions[action].signature.outputs,
                    f"{action} does not provide the required 2026.7 output")
        require(set(plugins["feature-table"].actions["summarize"].signature.outputs) ==
                {"feature_frequencies", "sample_frequencies", "summary"},
                "Unexpected feature-table summarize output contract")
        self.run("qiime-info", ["qiime", "info"])
        for binary in BIO_BINARIES[3:]:
            flag = "-version" if binary in ("blastn", "makeblastdb") else "--version"
            log = self.run("version-" + binary, [binary, flag])
            if binary in ("prefetch", "fasterq-dump"):
                require(sra.VERSION in log.read_text(), f"Unexpected {binary} version; expected {sra.VERSION}")
        for module in (None, "MetaDL", "AmpliconPIP", "AmpliconTAXA"):
            args = ["Meta2Data"] + ([module] if module else []) + ["--help"]
            self.run("help-" + (module or "Meta2Data"), args)
        self.save()

    def unit_tests(self):
        # Syntax checks do not create __pycache__ inside a read-only SIF.
        for path in [*sorted((ROOT / "scripts").glob("*.py")),
                     *sorted((ROOT / "tests").rglob("*.py"))]:
            compile(path.read_bytes(), str(path), "exec")
        for index, path in enumerate([*sorted((ROOT / "scripts").glob("*.sh")),
                                      *sorted((ROOT / "bin").glob("*"))]):
            self.run(f"shell-syntax-{index:02d}", ["bash", "-n", path])
        code = (
            "import json,sys,unittest\n"
            "suite=unittest.defaultTestLoader.discover('tests')\n"
            "result=unittest.TextTestRunner(verbosity=2).run(suite)\n"
            "summary=dict(tests=result.testsRun,failures=len(result.failures),"
            "errors=len(result.errors),skipped=result.skipped)\n"
            "open(sys.argv[1],'w').write(json.dumps(summary,indent=2)+'\\n')\n"
            "sys.exit(0 if result.wasSuccessful() and not result.skipped else 1)\n"
        )
        summary = self.work / "unit-tests.json"
        self.run("unit-tests", [sys.executable, "-c", code, summary])
        self.report["unit_tests"] = json.loads(summary.read_text())
        self.save()

    def prepare_references(self):
        # Fail before lengthy PIP runs if any release classifier is unavailable.
        sys.path.insert(0, str(ROOT / "scripts"))
        import reference_resources as refs
        cache_root = self.work / "references"
        kinds = ["gg2-sequences", "sepp", *(value[0] for value in CLASSIFIERS.values())]
        paths = {}
        for kind in kinds:
            self.run("reference-" + kind, [sys.executable, ROOT / "scripts/reference_resources.py",
                     "prepare", "--kind", kind, "--work-dir", cache_root])
            paths[kind] = cache_root / "results/db" / refs.RESOURCES[kind].filename
            require(paths[kind].is_file(), f"Reference was not published: {kind}")
        self.cache = cache_root / "results/db"
        self.report["references"] = {kind: str(path) for kind, path in paths.items()}
        self.save()
        return paths

    def case_directory(self, name):
        directory = self.work / name
        (directory / "results").mkdir(parents=True)
        # A shared fixture-only cache avoids copying large reference artifacts.
        (directory / "results/db").symlink_to(self.cache, target_is_directory=True)
        return directory

    def fixtures(self, reference):
        from generate_fixture import LOCAL_FIELDS, PACBIO_FORWARD, PACBIO_REVERSE, FORWARD, write_csv
        export = self.work / "reference-export"
        self.run("export-reference", ["qiime", "tools", "export", "--input-path", reference,
                                     "--output-path", export])
        fasta = export / "dna-sequences.fasta"
        require(fasta.is_file(), f"Missing exported GG2 sequences: {fasta}")
        fixture = self.work / "fixtures"
        self.run("generate-fixtures", [sys.executable, ROOT / "tests/integration/generate_fixture.py",
                 "--reference-fasta", fasta, "--output", fixture,
                 "--reads", 4000 if self.mode == "smoke" else 8000,
                 "--pacbio-reads", 1000 if self.mode == "smoke" else 2000,
                 "--templates", 3 if self.mode == "smoke" else 5])
        expected = json.loads((fixture / "expected_counts.json").read_text())["datasets"]
        rows = {}
        for filename in ("local_paired.csv", "local_single.csv", "local_pacbio.csv"):
            with (fixture / filename).open() as stream:
                rows.update((row["datasets"], row) for row in csv.DictReader(stream))
        if self.mode != "smoke":
            # These platform fixtures preserve the generated sequences/qualities.
            # They test executable/interface compatibility, not platform error models.
            for target, source, platform in (
                    ("ion", "se_explicit", "ION_TORRENT"),
                    ("ls454", "se_explicit", "LS454"),
                    ("ont", "pacbio", "OXFORD_NANOPORE")):
                target_dir = fixture / "local_data" / target
                target_dir.mkdir()
                samples = {}
                for index, (sample, record) in enumerate(sorted(expected[source]["samples"].items())):
                    new_sample = f"{target}_{chr(65 + index)}"
                    require(len(record["files"]) == 1, "Compatibility source must be single-end")
                    source_path = fixture / record["files"][0]["path"]
                    destination = target_dir / (new_sample + ".fastq.gz")
                    try:
                        os.link(source_path, destination)
                    except OSError:
                        shutil.copyfile(source_path, destination)
                    samples[new_sample] = {"expected_pipeline_reads": record["expected_pipeline_reads"]}
                expected[target] = {"platform": platform, "samples": samples}
                rows[target] = dict(datasets=target, path="local_data/" + target, platform=platform,
                                    primer_f=PACBIO_FORWARD if target == "ont" else FORWARD,
                                    primer_r=PACBIO_REVERSE if target == "ont" else "")
        names = ["pe_explicit"] if self.mode == "smoke" else list(rows)
        metadata = {}
        for method in ("dada2", "vsearch"):
            selected = [name for name in names if method == "vsearch" or name not in ("ls454", "ont")]
            path = fixture / f"compatibility-{method}.csv"
            write_csv(path, LOCAL_FIELDS, [rows[name] for name in selected])
            metadata[method] = (path, selected)
        return fixture, expected, metadata

    def load(self, path, semantic_type=None):
        from qiime2 import Artifact
        require(path.is_file(), f"Required artifact is missing: {path}")
        artifact = Artifact.load(path)
        artifact.validate(level="min")
        if semantic_type:
            require(str(artifact.type) == semantic_type,
                    f"{path.name}: expected {semantic_type}, got {artifact.type}")
        return artifact

    def table_and_sequences(self, table_path, sequences_path, expected_samples=None):
        from biom import Table
        import pandas as pd
        table = self.load(table_path, "FeatureTable[Frequency]").view(Table)
        sequences = self.load(sequences_path, "FeatureData[Sequence]").view(pd.Series)
        require(not sequences.index.has_duplicates, "Duplicate representative-sequence IDs")
        check_identifiers(table.ids(axis="observation"), sequences.index)
        values = table.matrix_data.data
        require(all(math.isfinite(float(v)) and v >= 0 and int(v) == v for v in values),
                f"Invalid abundance values in {table_path}")
        counts = dict(zip(map(str, table.ids(axis="sample")), map(int, table.sum(axis="sample"))))
        require(sum(counts.values()) > 0, f"No reads survived: {table_path}")
        if expected_samples is not None:
            require(set(counts) == set(expected_samples), f"Unexpected sample IDs: {table_path}: {counts}")
            for sample, count in counts.items():
                require(0 < count <= expected_samples[sample],
                        f"Read count for {sample} is outside 1..{expected_samples[sample]}: {count}")
        self.report["products"][str(table_path.relative_to(self.work))] = {
            "features": len(sequences), "sample_counts": counts}
        self.save()
        return table, sequences

    def pip(self, fixture, expected, metadata):
        from generate_fixture import fasta_records
        outputs = {}
        for method, (csv_path, names) in metadata.items():
            work = self.case_directory("pip-" + method)
            self.run("pip-" + method, ["Meta2Data", "AmpliconPIP",
                     "--local-m", csv_path, "--local-datasets-colNAME", "datasets",
                     "--local-path-colNAME", "path", "--local-platform-colNAME", "platform",
                     "--local-primer-f-colNAME", "primer_f", "--local-primer-r-colNAME", "primer_r",
                     "--" + method, "-t", self.threads], cwd=work)
            output = work / "results/pip"
            for name in names:
                samples = {sample: record["expected_pipeline_reads"]
                           for sample, record in expected[name]["samples"].items()}
                dataset = output / name
                _, sequences = self.table_and_sequences(
                    dataset / f"{name}-{method}-final-table.qza",
                    dataset / f"{name}-{method}-final-rep-seqs.qza", samples)
                if name == "pacbio" and method == "dada2":
                    inserts = {seq for _, seq in fasta_records(fixture / "pacbio_inserts.fasta")}
                    require(set(map(str, sequences)) <= inserts,
                            "PacBio DADA2 representatives do not match the known primer-free inserts")
            outputs[method] = output
        return outputs

    def summary(self, path):
        from qiime2 import Visualization
        require(path.is_file(), f"Feature-table summary was not produced: {path}")
        Visualization.load(path).validate(level="min")
        for suffix in ("feature-frequencies", "sample-frequencies"):
            self.load(path.with_name(path.stem + "-" + suffix + ".qza"), "ImmutableMetadata")

    def taxa(self, inputs):
        import pandas as pd
        from skbio import TreeNode
        sys.path.insert(0, str(ROOT / "scripts"))
        import reference_resources as refs
        # Reuse the small PE VSEARCH output for two additional tree branches.
        cases = [(key, "notree", inputs["dada2"], "dada2") for key in CLASSIFIERS]
        cases += [("gg2-2024.09", tree, inputs["vsearch"] / "pe_explicit", "vsearch")
                  for tree in ("sepp", "singleV")]
        for classifier, tree, source, method in cases:
            name = "taxa-" + classifier + "-" + tree
            work = self.case_directory(name)
            command = ["Meta2Data", "AmpliconTAXA", "-i", source, "--" + method,
                       "--classifier", classifier, "-t", self.threads]
            if tree != "sepp":
                command += ["--" + tree]
            self.run(name, command, cwd=work)
            label = CLASSIFIERS[classifier][1]
            files = list((work / "results/taxa").rglob(label + "Taxonomy.qza"))
            require(len(files) == 1, f"Expected one taxonomy artifact in {work}, found {files}")
            taxonomy_path = files[0]
            final = taxonomy_path.parent
            taxonomy = self.load(taxonomy_path, "FeatureData[Taxonomy]").view(pd.DataFrame)
            require({"Taxon", "Confidence"} <= set(taxonomy), "Taxonomy columns changed")
            require(all(math.isfinite(float(v)) and 0 <= float(v) <= 1 for v in taxonomy["Confidence"]),
                    "Invalid taxonomy confidence values")
            resource = refs.RESOURCES[CLASSIFIERS[classifier][0]]
            identity = check_taxonomy_classifier(taxonomy_path, self.cache / resource.filename)
            ranks = (check_silva144_ranks(taxonomy["Taxon"]) if classifier == "silva-144" else
                     {"observed_rank_counts": sorted(
                         {len(str(taxon).split(";")) for taxon in taxonomy["Taxon"]})})
            if tree == "sepp":
                table_path, sequences_path = final / "treeFilteredTable.qza", final / "treeFilteredRepSeqs.qza"
            else:
                table_path, sequences_path = final / "orientedTable.qza", final / "orientedRepSeqs.qza"
            table, sequences = self.table_and_sequences(table_path, sequences_path)
            check_identifiers(table.ids(axis="observation"), sequences.index, taxonomy.index)
            require(not taxonomy.index.has_duplicates, "Duplicate feature IDs in taxonomy")
            if tree != "sepp":
                require(set(taxonomy.index) == set(sequences.index), "Taxonomy and oriented feature IDs differ")
            if tree != "notree":
                tree_path = final / ("seppTree.qza" if tree == "sepp" else "denovoRootedTree.qza")
                tips = {tip.name for tip in self.load(tree_path, "Phylogeny[Rooted]").view(TreeNode).tips()}
                check_identifiers(table.ids(axis="observation"), sequences.index, taxonomy.index, tips)
                if tree == "singleV":
                    require(tips == set(sequences.index), "De novo tree has unexpected tips")
            self.summary(final / "tmp/mergedTableSummary.qzv")
            self.summary(final / "tmp" / ("treeFilteredTableSummary.qzv" if tree == "sepp"
                                          else "orientedTableSummary.qzv"))
            self.report["products"][str(taxonomy_path.relative_to(self.work))] = {
                "features": len(taxonomy), "requested_classifier": classifier,
                **identity, **ranks,
            }
            self.save()

    def execute(self):
        self.environment()
        if self.mode in ("full", "checks"):
            self.unit_tests()
        if self.mode in ("full", "smoke"):
            references = self.prepare_references()
            fixture, expected, metadata = self.fixtures(references["gg2-sequences"])
            self.taxa(self.pip(fixture, expected, metadata))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True,
                        help="New or empty directory for evidence, fixtures and results")
    parser.add_argument("--threads", type=int, default=2)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--checks-only", action="store_true", help="Environment and complete unit tests only")
    modes.add_argument("--environment-only", action="store_true", help="Runtime interfaces and entrypoints only")
    modes.add_argument("--smoke", action="store_true", help="Packaged-image fixtures; unit/platform matrix already validated")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    mode = ("environment" if args.environment_only else "checks" if args.checks_only
            else "smoke" if args.smoke else "full")
    validation = None
    try:
        validation = Validation(args.work_dir.resolve(), mode, args.threads)
        validation.execute()
    except Exception as error:
        if validation is not None:
            validation.report.update(status="failed", error=str(error),
                                     finished=datetime.now(timezone.utc).isoformat())
            validation.save()
        print(f"Validation failed: {error}", file=sys.stderr)
        return 1
    validation.report.update(status="passed", finished=datetime.now(timezone.utc).isoformat())
    validation.save()
    print(f"Validation passed ({mode}); evidence: {validation.work / 'validation.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
