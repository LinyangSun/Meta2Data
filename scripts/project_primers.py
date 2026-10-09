"""Validate project-specific primers and prepare a metadata subset before running PIP."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import sys
import tempfile


PROJECT = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
RUN = re.compile(r"(?:[EDS]RR|CRR)[0-9]+\Z")
PRIMER = re.compile(r"[ACGTRYSWKMBDHVN]+\Z", re.IGNORECASE)


def validate_primers(projects, forward, reverse=None):
    """Return an ordered mapping without broadcasting or filling missing primers."""
    reverse = reverse or []
    if not projects:
        raise ValueError("--public-bioprojectIDs requires at least one project")
    if any(not PROJECT.fullmatch(project) for project in projects):
        raise ValueError("Project IDs must contain only letters, digits, dots, underscores or hyphens")
    if len(set(projects)) != len(projects):
        raise ValueError("--public-bioprojectIDs contains duplicate project IDs")
    if len(forward) != len(projects):
        raise ValueError("--public-primer-fwd must provide one primer per --public-bioprojectIDs entry")
    if reverse and len(reverse) != len(projects):
        raise ValueError("--public-primer-rev must provide one primer per --public-bioprojectIDs entry, or be omitted entirely")
    for value in [*forward, *reverse]:
        if not PRIMER.fullmatch(value):
            raise ValueError(f"Invalid primer sequence {value!r}; use IUPAC DNA bases, not names or paths")
    return {project: {"forward": forward[index].upper(),
                      "reverse": reverse[index].upper() if reverse else ""}
            for index, project in enumerate(projects)}


def _atomic_write(destination, write):
    destination = Path(destination)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                     prefix=f".{destination.name}.", dir=destination.parent,
                                     delete=False) as handle:
        temporary = Path(handle.name)
        try:
            write(handle)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def prepare(metadata, output, col_bioproject, col_sra, projects, forward, reverse=None):
    mapping = validate_primers(projects, forward, reverse)
    source = Path(metadata).resolve()
    destination = Path(output).resolve()
    selected_path = destination / "selected_metadata.csv"
    map_path = destination / "project_primers.json"
    if source in (selected_path, map_path):
        raise ValueError("Input metadata must differ from the generated selected_metadata.csv; place the input metadata outside results/pip")
    selected = []
    seen = set()
    counts = {project: 0 for project in projects}
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        for column in (col_bioproject, col_sra):
            if column not in fields:
                raise ValueError(f"Metadata is missing column {column!r}")
        if col_bioproject == col_sra:
            raise ValueError("Project and Run column names must differ")
        for row in reader:
            project = (row.get(col_bioproject) or "").strip()
            if project not in mapping:
                continue
            seen.add(project)
            run = (row.get(col_sra) or "").strip()
            if not run or run.lower() in {"na", "nan", "none", "null"}:
                continue
            if not RUN.fullmatch(run):
                raise ValueError(f"Invalid Run {run!r} in project {project}; expected SRR, ERR, DRR or CRR accession")
            row[col_bioproject], row[col_sra] = project, run
            selected.append(row)
            counts[project] += 1
    missing = [project for project in projects if project not in seen]
    if missing:
        raise ValueError("Projects not present in metadata: " + ", ".join(missing))
    empty = [project for project in projects if not counts[project]]
    if empty:
        raise ValueError("Projects have no usable Run accessions: " + ", ".join(empty))
    record = {"schema_version": 1, "source_metadata": str(source),
              "project_column": col_bioproject, "run_column": col_sra,
              "project_order": projects, "projects": mapping}
    # All validation completes before replacing any output files.
    destination.mkdir(parents=True, exist_ok=True)

    def write_csv(handle):
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)

    _atomic_write(selected_path, write_csv)
    _atomic_write(map_path, lambda handle: handle.write(json.dumps(record, indent=2) + "\n"))
    return selected_path, map_path


def read_mapping(path):
    with Path(path).open(encoding="utf-8") as handle:
        record = json.load(handle)
    if record.get("schema_version") != 1:
        raise ValueError("Unsupported project-primer mapping version")
    order = record["project_order"]
    projects = record["projects"]
    if not isinstance(order, list) or not isinstance(projects, dict) or set(order) != set(projects):
        raise ValueError("Invalid project-primer mapping")
    forwards = [projects[project]["forward"] for project in order]
    reverses = [projects[project]["reverse"] for project in order]
    # Empty reverse lists mean explicit forward-only processing for every project.
    if all(not reverse for reverse in reverses):
        reverses = []
    return validate_primers(order, forwards, reverses)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--metadata", required=True)
    prepare_parser.add_argument("--output", required=True)
    prepare_parser.add_argument("--col-bioproject", required=True)
    prepare_parser.add_argument("--col-sra", required=True)
    prepare_parser.add_argument("--projects", nargs="+", required=True)
    prepare_parser.add_argument("--forward", nargs="+", required=True)
    prepare_parser.add_argument("--reverse", nargs="+", default=[])
    emit_parser = commands.add_parser("emit")
    emit_parser.add_argument("--mapping", required=True)
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            prepare(args.metadata, args.output, args.col_bioproject, args.col_sra,
                    args.projects, args.forward, args.reverse)
        else:
            for project, primers in read_mapping(args.mapping).items():
                print(project, primers["forward"], primers["reverse"], sep="\t")
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()
