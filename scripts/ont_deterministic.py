#!/usr/bin/env python3
"""Internal, lossless ordering of ONT intermediates before order-sensitive tools.

GNU sort provides a disk-backed, bytewise total order with a 64 MiB sort buffer.
Only one record is parsed at a time. Sequence case and orientation, duplicate
records, FASTQ quality, and complete FASTA headers (including size annotations)
are preserved. FASTQ-to-FASTA conversion deliberately replaces original names:
sample hash + sequence hash + sequence rank + duplicate ordinal is unique
even for duplicate input names or a sequence-hash collision.

Sorting permits input == output: the input is fully read and closed before an
atomic replacement. Conversion to FASTA refuses to overwrite its FASTQ input.
Gzip inputs are accepted; output paths must be uncompressed. All temporary
files are placed beside the output and removed on failure.
"""

import argparse
import base64
import gzip
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile


DNA = frozenset("ACGTURYKMSWBDHVNacgturykmswbdhvn")


def _open_input(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="ascii", newline=None)


def _line(stream):
    value = stream.readline()
    return None if value == "" else value.rstrip("\r\n")


def _sequence(value, label):
    if not value or not set(value) <= DNA:
        raise ValueError(f"{label}: empty sequence or invalid nucleotide characters")
    return value


def _fastq_records(stream):
    """Read complete FASTQ records, including wrapped sequence/quality lines."""
    count = 0
    while True:
        header = _line(stream)
        if header is None:
            return
        count += 1
        label = f"FASTQ record {count}"
        if not header.startswith("@") or not header[1:].split():
            raise ValueError(f"{label}: missing read name")
        parts = []
        while True:
            line = _line(stream)
            if line is None:
                raise ValueError(f"{label}: truncated sequence or missing '+' line")
            if line.startswith("+"):
                plus = line
                break
            parts.append(_sequence(line, label))
        sequence = _sequence("".join(parts), label)
        if plus[1:].strip() and plus[1:].split()[0] != header[1:].split()[0]:
            raise ValueError(f"{label}: '+' identifier differs from read name")
        quality_parts, length = [], 0
        while length < len(sequence):
            line = _line(stream)
            if line is None or not line:
                raise ValueError(f"{label}: truncated or empty quality line")
            if any(ord(char) < 33 or ord(char) > 126 for char in line):
                raise ValueError(f"{label}: invalid quality characters")
            quality_parts.append(line)
            length += len(line)
        if length != len(sequence):
            raise ValueError(f"{label}: sequence and quality lengths differ")
        yield sequence, "".join(quality_parts), header, plus


def _fasta_records(stream):
    header, parts = None, []
    for raw in stream:
        line = raw.rstrip("\r\n")
        if line.startswith(">"):
            if header is not None:
                yield _sequence("".join(parts), "FASTA record"), header
            if not line[1:].split():
                raise ValueError("FASTA record: missing sequence name")
            header, parts = line, []
        elif header is None:
            raise ValueError("FASTA data before first header")
        else:
            parts.append(_sequence(line, "FASTA record"))
    if header is not None:
        yield _sequence("".join(parts), "FASTA record"), header


def _encode(value):
    return base64.b64encode(value.encode("ascii")).decode("ascii")


def _decode(value):
    return base64.b64decode(value, validate=True).decode("ascii")


def _same_file(source, destination):
    return (Path(source).resolve() == Path(destination).resolve()
            or (Path(destination).exists() and os.path.samefile(source, destination)))


def _transform(source, destination, mode, sample=None):
    source, destination = Path(source), Path(destination)
    if destination.suffix == ".gz":
        raise ValueError("Ordered intermediates require an uncompressed output path")
    if mode == "to-fasta":
        if not sample:
            raise ValueError("A nonempty sample namespace is required")
        if _same_file(source, destination):
            raise ValueError("FASTQ-to-FASTA conversion cannot overwrite its input")
        # A bounded namespace keeps even long sample names within SAM QNAME
        # limits. Domain separation distinguishes this hash from sequence IDs.
        namespace = hashlib.sha256(b"meta2data-ont-sample\0" + sample.encode("utf-8")).hexdigest()
    # Intermediate files must be accessible through the same container bind as
    # the output, and os.replace must stay on the same filesystem.
    with tempfile.TemporaryDirectory(prefix=".ont-sort-", dir=destination.parent) as folder:
        work = Path(folder)
        unsorted, ordered, result = (work / name for name in ("records", "ordered", "result"))
        with _open_input(source) as stream, unsorted.open("w", encoding="ascii", newline="\n") as sink:
            if mode == "sort-fasta":
                for sequence, header in _fasta_records(stream):
                    sink.write(f"{sequence}\t{_encode(header)}\n")
            else:
                for sequence, quality, header, plus in _fastq_records(stream):
                    if mode == "to-fasta":
                        sink.write(sequence + "\n")
                    else:
                        sink.write(f"{sequence}\t{quality}\t{_encode(header)}\t{_encode(plus)}\n")
        with ordered.open("wb") as sink:
            subprocess.run(
                ["sort", "--buffer-size=64M", "--parallel=1",
                 "--temporary-directory", str(work), "--", str(unsorted)],
                stdout=sink, check=True, env={**os.environ, "LC_ALL": "C"},
            )
        with ordered.open("r", encoding="ascii") as stream, result.open("w", encoding="ascii", newline="\n") as sink:
            previous, rank, duplicate = None, 0, 0
            for raw in stream:
                fields = raw.rstrip("\n").split("\t")
                sequence = fields[0]
                if mode == "to-fasta":
                    if sequence != previous:
                        rank += 1
                        duplicate = 0
                        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
                        previous = sequence
                    duplicate += 1
                    sink.write(f">s{namespace}.{digest}.{rank}.{duplicate}\n{sequence}\n")
                elif mode == "sort-fasta":
                    sink.write(f"{_decode(fields[1])}\n{sequence}\n")
                else:
                    sink.write(f"{_decode(fields[2])}\n{sequence}\n{_decode(fields[3])}\n{fields[1]}\n")
        # No output is changed until input validation and external sorting have
        # both succeeded. Closing result above also detects buffered write errors.
        os.replace(result, destination)


def sort_fastq(source, destination):
    _transform(source, destination, "sort-fastq")


def fastq_to_fasta(source, destination, sample):
    _transform(source, destination, "to-fasta", sample)


def sort_fasta(source, destination):
    _transform(source, destination, "sort-fasta")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("sort-fastq", "to-fasta", "sort-fasta"):
        command = commands.add_parser(name)
        command.add_argument("--input", required=True)
        command.add_argument("--output", required=True)
        if name == "to-fasta":
            command.add_argument("--sample", required=True)
    args = parser.parse_args()
    try:
        _transform(args.input, args.output, args.command, getattr(args, "sample", None))
    except (OSError, EOFError, UnicodeError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"[ERROR] ONT deterministic ordering failed: {error}\n")


if __name__ == "__main__":
    main()
