"""Synchronize R2's 5' primer cut only after R1's cut is authorized.

Automatic callers must run their existing invalid/unknown-primer gates before
calling synchronize_plan. Explicit primers are detected by cutadapt without
changing reads; the effective cuts are then applied once to original mates.
"""
from __future__ import annotations

from collections import Counter
import os
from pathlib import Path
import subprocess
import tempfile

from illumina_pe import open_fastq, paired_records, read_id, write_record


def synchronize_plan(info, forward_trim, reverse_trim, trims=None):
    """Keep independent detection evidence while overriding effective R2 cuts."""
    if info.get('detect_only') or info.get('layout') != 'PE':
        return reverse_trim
    info['r2_sync_policy'] = 'r1_authorized_5prime_trim'
    info['r2_synchronized'] = False
    if trims is not None:
        for r1_group, r2_group in (('forward', 'reverse'), ('reverse', 'forward')):
            length = trims['r1'][r1_group]
            if length > 0:
                entry = info['original_ends']['r2'][r2_group]
                entry['independent_trim_length'] = trims['r2'][r2_group]
                entry['trim_source'] = 'R1/' + r1_group
                trims['r2'][r2_group] = length
                info['r2_synchronized'] = True
    elif forward_trim > 0:
        entry = info['reverse_primer']
        entry['independent_trim_length'] = reverse_trim
        entry['trim_source'] = 'R1'
        reverse_trim = forward_trim
        info['r2_synchronized'] = True
    return reverse_trim


def info_cut(line, original):
    """Read a cutadapt info row for one regular 5' adapter (default one match)."""
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 2 or read_id(fields[0]) != read_id(original[0]):
        raise ValueError('cutadapt info read IDs/order do not match the original FASTQ')
    if fields[1] == '-1':
        return 0
    if len(fields) < 8:
        raise ValueError('Incomplete cutadapt info row')
    start, end = int(fields[2]), int(fields[3])
    if not 0 <= start < end <= len(original[1]):
        raise ValueError('Invalid cutadapt 5-prime match coordinates')
    return end


def explicit_paired(r1, r2, output1, output2, forward, reverse='', runner=subprocess.run):
    """Match each supplied primer, decide both cuts, then write mates once."""
    if not forward:
        raise ValueError('A forward primer is required for explicit paired synchronization')
    out1, out2 = Path(output1), Path(output2)
    if out1.resolve() in (Path(r1).resolve(), Path(r2).resolve()) or out2.resolve() in (
            Path(r1).resolve(), Path(r2).resolve()) or out1.resolve() == out2.resolve():
        raise ValueError('Primer outputs must be distinct from the original mates')
    result = dict(pairs=0, synchronized_pairs=0, r1_trim_lengths=Counter(),
                  r2_independent_trim_lengths=Counter(), r2_effective_trim_lengths=Counter(),
                  commands=[], detection_reports=[])
    with tempfile.TemporaryDirectory(prefix='.primer-sync-', dir=out1.parent) as temporary:
        temporary = Path(temporary)
        info_paths = []
        for index, source, primer in ((1, r1, forward), (2, r2, reverse)):
            if not primer:
                info_paths.append(None)
                continue
            path = temporary / f'R{index}.info.tsv'
            log = temporary / f'R{index}.cutadapt.log'
            command = ['cutadapt', '--action=none', '--info-file', str(path),
                       '-g', primer, '-o', os.devnull, str(source)]
            result['commands'].append(command)
            with log.open('w') as stream:
                runner(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
            result['detection_reports'].append(log.read_text(errors='replace'))
            info_paths.append(path)
        # ExitStack closes every stream if a mate or detector output is corrupt.
        from contextlib import ExitStack
        with ExitStack() as stack:
            first_info = stack.enter_context(info_paths[0].open())
            second_info = stack.enter_context(info_paths[1].open()) if info_paths[1] else None
            first_out = stack.enter_context(open_fastq(out1, 'wt'))
            second_out = stack.enter_context(open_fastq(out2, 'wt'))
            for original1, original2 in paired_records(r1, r2):
                first_line = first_info.readline()
                second_line = second_info.readline() if second_info else None
                if not first_line or (second_info is not None and not second_line):
                    raise ValueError('cutadapt info contains fewer reads than the paired input')
                trim1 = info_cut(first_line, original1)
                independent2 = info_cut(second_line, original2) if second_line is not None else 0
                trim2 = trim1 if trim1 > 0 else independent2
                result['pairs'] += 1
                result['synchronized_pairs'] += trim1 > 0
                result['r1_trim_lengths'][trim1] += 1
                result['r2_independent_trim_lengths'][independent2] += 1
                result['r2_effective_trim_lengths'][trim2] += 1
                for stream, original, length in ((first_out, original1, trim1),
                                                  (second_out, original2, trim2)):
                    header, sequence, plus, quality = original
                    write_record(stream, (header, sequence[length:], plus, quality[length:]))
            if first_info.readline() or (second_info is not None and second_info.readline()):
                raise ValueError('cutadapt info contains more reads than the paired input')
    return result
