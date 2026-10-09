#!/usr/bin/env python3
"""Generate reproducible synthetic amplicon reads for offline integration checks.

Templates come from a supplied GG2 FASTA. Accessions in the generated public
metadata are fictional: use public_fixture.json with the archive test shim,
never query real archives for these identifiers. Sequencing errors are sampled
from the Phred probabilities; these fixtures do not establish biological accuracy.
"""
import argparse
from collections import Counter
import csv
import difflib
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import random
import re
import shutil

FORWARD = 'CCTACGGGNGGCWGCAG'
REVERSE = 'GACTACHVGGGTATCTAATCC'
PACBIO_FORWARD = 'AGAGTTTGATCCTGGCTCAG'
PACBIO_REVERSE = 'TACGGYTACCTTGTTAYGACTT'
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT',
         'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
         'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
ALTERNATIVES = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
LOCAL_FIELDS = ['datasets', 'path', 'platform', 'primer_f', 'primer_r']
PROJECT = 'PRJNA999990001'
RUNS = ['SRR999990001', 'SRR999990002']


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def primer_pattern(primer):
    return re.compile(''.join(base if len(IUPAC[base]) == 1 else '[' + IUPAC[base] + ']'
                              for base in primer))


def fasta_records(path):
    name, chunks = None, []
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks).upper()


def build_pacbio_template(full):
    """Replace the native terminal 1492R binding site, never duplicate it.

    Some GG2 records stop after the first 19 bases of the reverse-complemented
    primer. Their missing terminal bases are supplied by the experimental
    primer. A complete native site may have downstream sequence, which is
    outside the PCR product and must not remain in the synthetic amplicon.
    """
    reverse_target = reverse_complement(PACBIO_REVERSE)
    core_length = 19
    core = primer_pattern(reverse_target[:core_length])
    sites = list(core.finditer(full))
    if len(sites) != 1:
        raise ValueError('1492R requires one unique native binding site')
    site = sites[0]
    if site.start() < max(len(PACBIO_FORWARD), len(full) - 100):
        raise ValueError('1492R binding site is outside the terminal 100 bp')
    native_end = min(site.start() + len(reverse_target), len(full))
    covered = full[site.start():native_end]
    if not primer_pattern(reverse_target[:len(covered)]).fullmatch(covered):
        raise ValueError('1492R native binding sequence is incompatible with the primer')
    insert = full[len(PACBIO_FORWARD):site.start()]
    concrete_reverse = ''.join(IUPAC[base][0] for base in PACBIO_REVERSE)
    pacbio = PACBIO_FORWARD + insert + reverse_complement(concrete_reverse)
    if len(pacbio) <= 1400:
        raise ValueError('1492R-defined PacBio amplicon does not exceed 1400 bp')
    # This also rejects artificial repeated sites in the input references.
    if len(list(core.finditer(pacbio))) != 1:
        raise ValueError('1492R construction produced multiple binding sites')
    binding = dict(coordinate_system='0-based, half-open reference coordinates',
                   forward_replaced_start=0, forward_replaced_end=len(PACBIO_FORWARD),
                   reverse_start=site.start(), reverse_covered_end=native_end,
                   reverse_reference_bases=len(covered), reverse_reference_sequence=covered,
                   downstream_reference_bases_removed=len(full) - native_end,
                   insert_start=len(PACBIO_FORWARD), insert_end=site.start())
    return dict(pacbio=pacbio, pacbio_insert=insert, pacbio_binding=binding)


def select_templates(reference, count):
    forward = primer_pattern(FORWARD)
    reverse = primer_pattern(reverse_complement(REVERSE))
    selected = []
    pacbio_rejections = Counter()
    for name, full in fasta_records(reference):
        # Leave enough sequence after CCS primer removal to remain full length.
        if not 1440 <= len(full) <= 1510 or set(full) - set('ACGT'):
            continue
        left = forward.search(full)
        right = reverse.search(full, left.end()) if left else None
        if not left or not right or not 250 <= left.start() <= 400:
            continue
        amplicon = full[left.start():right.end()]
        if not 400 <= len(amplicon) <= 510:
            continue
        # Retain distinct templates so a 97% clustering run has several features.
        if any(difflib.SequenceMatcher(None, amplicon, previous['amplicon'],
                                      autojunk=False).ratio() > 0.94
               for previous in selected):
            continue
        try:
            pacbio = build_pacbio_template(full)
        except ValueError as error:
            pacbio_rejections[str(error)] += 1
            continue
        selected.append(dict(id=f'template_{len(selected) + 1}', reference_id=name,
                             full=full, amplicon=amplicon, **pacbio,
                             forward_start=left.start(), reverse_end=right.end()))
        if len(selected) == count:
            return selected
    raise ValueError(f'Found only {len(selected)} sufficiently distinct templates; need {count}; '
                     f'PacBio binding-site rejections: {dict(pacbio_rejections)}')


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def write_csv(path, fields, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def gzip_writer(path):
    # Blank gzip filename and mtime=0 make bytes/checksums deterministic.
    return io.TextIOWrapper(gzip.GzipFile(filename='', mode='wb', compresslevel=6,
                                        fileobj=Path(path).open('wb'), mtime=0),
                            encoding='ascii', newline='\n')


def digest(path, algorithm='md5'):
    result = hashlib.new(algorithm)
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def rng_for(seed, label):
    key = hashlib.sha256(f'{seed}:{label}'.encode()).digest()
    return random.Random(int.from_bytes(key, 'big'))


def observed_read(sequence, rng, quality_counts, error_counts, protect_both=False):
    bases, qualities = [], []
    last = max(len(sequence) - 1, 1)
    for index, base in enumerate(sequence):
        quality = max(28, min(40, round(38 - 8 * index / last) + rng.randint(-3, 3)))
        quality_counts[str(quality)] += 1
        qualities.append(chr(quality + 33))
        protected = index < 24 if protect_both else index < 20
        protected = protected or (protect_both and index >= len(sequence) - 24)
        if not protected and rng.random() < 10 ** (-quality / 10):
            base = ALTERNATIVES[base][rng.randrange(3)]
            error_counts[str(quality)] += 1
        bases.append(base)
    return ''.join(bases), ''.join(qualities)


def fastq_record(name, sequence, rng, quality_counts, error_counts, protect_both=False):
    read, quality = observed_read(sequence, rng, quality_counts, error_counts, protect_both)
    return f'@{name}\n{read}\n+\n{quality}\n'


def generate_sample(output, dataset, sample, templates, count, seed, paired=False,
                    chunks=False, pacbio=False, sample_index=0):
    folder = output / 'local_data' / dataset
    folder.mkdir(parents=True, exist_ok=True)
    rng = rng_for(seed, dataset + '/' + sample)
    weights = [40, 25, 18, 10, 5, 2][:len(templates)]
    weights = weights[sample_index:] + weights[:sample_index]
    selected = rng.choices(range(len(templates)), weights=weights, k=count)
    template_counts = Counter(templates[index]['id'] for index in selected)
    quality_counts, error_counts = Counter(), Counter()
    file_info = []
    shard_count = min(4, count) if chunks else 1
    start = 0
    for shard in range(shard_count):
        size = count // shard_count + (shard < count % shard_count)
        suffix = f'_L{shard // 2 + 1:03d}_R{{mate}}_{shard % 2 + 1:03d}' if chunks else '_R{mate}'
        paths = [folder / (sample + suffix.format(mate=mate) + '.fastq.gz')
                 for mate in (1, 2)] if paired else [folder / (sample + '.fastq.gz')]
        streams = [gzip_writer(path) for path in paths]
        try:
            for index in range(start, start + size):
                template = templates[selected[index]]
                if pacbio:
                    seq = template['pacbio']
                    if index % 2:
                        seq = reverse_complement(seq)
                    sequences = [seq]
                elif paired:
                    amplicon = template['amplicon']
                    sequences = [amplicon[:300], reverse_complement(amplicon[-300:])]
                else:
                    sequences = [template['amplicon'][:260]]
                for mate, (stream, sequence) in enumerate(zip(streams, sequences), 1):
                    # Match pair IDs after stripping the optional /1 and /2 suffix.
                    name = f'{sample}.{index + 1}' + (f'/{mate}' if paired else '')
                    stream.write(fastq_record(name, sequence, rng, quality_counts,
                                              error_counts, protect_both=pacbio))
        finally:
            for stream in streams:
                stream.close()
        for mate, path in enumerate(paths, 1):
            file_info.append(dict(path=path.relative_to(output).as_posix(), mate=mate,
                                  records=size, md5=digest(path), bytes=path.stat().st_size))
        start += size
    return dict(sample_id=sample, layout='PAIRED' if paired else 'SINGLE',
                expected_pipeline_reads=count, individual_read_records=count * (2 if paired else 1),
                read_length=(None if pacbio else 300 if paired else 260),
                template_counts=dict(template_counts), files=file_info,
                phred_base_counts=dict(quality_counts), injected_errors_by_phred=dict(error_counts))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reads', type=int, default=8000,
                        help='Pairs or single-end reads per sample (default: 8000)')
    parser.add_argument('--pacbio-reads', type=int, default=2000)
    parser.add_argument('--templates', type=int, default=5, choices=range(3, 7))
    parser.add_argument('--seed', type=int, default=20261009)
    args = parser.parse_args()
    if args.reads < 4 or args.pacbio_reads < 1:
        parser.error('--reads must be at least 4 and --pacbio-reads positive')
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('--output must be absent or empty; existing fixtures are never overwritten')
    templates = select_templates(args.reference_fasta, args.templates)
    output.mkdir(parents=True, exist_ok=True)
    for key, filename in [('full', 'reference_templates.fasta'),
                          ('amplicon', 'amplicon_templates.fasta'),
                          ('pacbio', 'pacbio_templates.fasta'),
                          ('pacbio_insert', 'pacbio_inserts.fasta')]:
        (output / filename).write_text(''.join(f">{item['id']} {item['reference_id']}\n{item[key]}\n"
                                             for item in templates))
    expected = dict(schema_version=1, synthetic=True, seed=args.seed, datasets={})
    specifications = [('pe_auto', 'autoPE', True, True, False),
                      ('pe_explicit', 'explicitPE', True, False, False),
                      ('se_auto', 'autoSE', False, False, False),
                      ('se_explicit', 'explicitSE', False, False, False),
                      ('pacbio', 'ccs', False, False, True)]
    for dataset, prefix, paired, chunks, pacbio in specifications:
        samples = {}
        for index, suffix in enumerate(('A', 'B')):
            sample = prefix + '_' + suffix
            samples[sample] = generate_sample(output, dataset, sample, templates,
                                             args.pacbio_reads if pacbio else args.reads,
                                             args.seed, paired, chunks, pacbio, index)
        expected['datasets'][dataset] = dict(platform='PACBIO_SMRT' if pacbio else 'ILLUMINA',
            samples=samples, expected_pipeline_reads=sum(s['expected_pipeline_reads'] for s in samples.values()))
        print(f'Generated {dataset}: {len(samples)} samples', flush=True)

    local_rows = {}
    for dataset, _, _, _, pacbio in specifications:
        explicit = dataset.endswith('_explicit') or pacbio
        local_rows[dataset] = dict(datasets=dataset, path='local_data/' + dataset,
            platform='PACBIO_SMRT' if pacbio else 'ILLUMINA',
            primer_f=PACBIO_FORWARD if pacbio else FORWARD if explicit else '',
            primer_r=PACBIO_REVERSE if pacbio else REVERSE if dataset == 'pe_explicit' else '')
    for filename, names in [('local_all.csv', ['pe_auto', 'se_explicit']),
                            ('local_paired.csv', ['pe_auto', 'pe_explicit']),
                            ('local_single.csv', ['se_auto', 'se_explicit']),
                            ('local_pacbio.csv', ['pacbio'])]:
        write_csv(output / filename, LOCAL_FIELDS, [local_rows[name] for name in names])

    archive = output / 'archive'
    archive.mkdir()
    public_runs = {}
    for run, sample in zip(RUNS, ('explicitPE_A', 'explicitPE_B')):
        record = expected['datasets']['pe_explicit']['samples'][sample]
        files = []
        for item in record['files']:
            destination = archive / f"{run}_{item['mate']}.fastq.gz"
            try:
                os.link(output / item['path'], destination)
            except OSError:
                shutil.copyfile(output / item['path'], destination)
            files.append(dict(path=destination.relative_to(output).as_posix(),
                              url='ftp.sra.ebi.ac.uk/fixtures/' + destination.name,
                              md5=digest(destination), bytes=destination.stat().st_size,
                              mate=item['mate'], records=item['records']))
        public_runs[run] = dict(platform='ILLUMINA', layout='PAIRED', read_pairs=args.reads,
                               files=files, template_counts=record['template_counts'])
    write_csv(output / 'online_metadata.csv', ['Bioproject', 'Run'],
              [dict(Bioproject=PROJECT, Run=run) for run in RUNS])
    manifest = dict(schema_version=1, synthetic=True, seed=args.seed,
        warning='Fictional accessions and archive URLs: use the integration shim; do not query real archives.',
        metadata='online_metadata.csv', projects={PROJECT: dict(platform='ILLUMINA', runs=RUNS)},
        runs=public_runs, primers=dict(forward=FORWARD, reverse=REVERSE,
                                     pacbio_forward=PACBIO_FORWARD, pacbio_reverse=PACBIO_REVERSE),
        reference_fasta=str(args.reference_fasta.resolve()),
        templates=[dict(id=t['id'], reference_id=t['reference_id'], full_length=len(t['full']),
                        amplicon_length=len(t['amplicon']), pacbio_length=len(t['pacbio']),
                        forward_start=t['forward_start'], reverse_end=t['reverse_end'],
                        pacbio_binding=t['pacbio_binding'],
                        pacbio_insert_length=len(t['pacbio_insert']),
                        pacbio_insert_sha256=hashlib.sha256(t['pacbio_insert'].encode()).hexdigest(),
                        full_sha256=hashlib.sha256(t['full'].encode()).hexdigest()) for t in templates],
        template_files=dict(reference='reference_templates.fasta', amplicon='amplicon_templates.fasta',
                            pacbio='pacbio_templates.fasta', pacbio_insert='pacbio_inserts.fasta'),
        expected_counts='expected_counts.json',
        simulation=dict(quality_range=[28, 40], substitution_probability='10 ** (-Q / 10)',
                        protected_prefix_bases=20, paired_read_length=300, single_read_length=260,
                        pacbio_protected_end_bases=24, pacbio_reverse_orientation_fraction=0.5,
                        pacbio_construction='Replace leading reference bases with explicit 27F; retain the insert before the unique terminal native 1492R binding site; replace that site with the reverse complement of concrete 1492R and remove downstream reference sequence.'))
    write_json(output / 'public_fixture.json', manifest)
    write_json(output / 'expected_counts.json', expected)
    (output / 'SYNTHETIC_DATA.txt').write_text(
        'Synthetic integration fixtures; not an experimental dataset.\n'
        'Archive IDs and URLs are fictional and require the integration archive shim.\n'
        'FASTQ sequences derive from the listed GG2 references with seeded substitutions.\n'
        'Identical seed, templates and read counts produce identical compressed FASTQ bytes.\n'
        'File sample IDs are unique across local datasets; archive files reuse explicit PE content.\n')
    print(output / 'public_fixture.json', flush=True)


if __name__ == '__main__':
    main()
