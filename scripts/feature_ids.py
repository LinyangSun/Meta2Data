"""Give vsearch features stable sequence IDs before merging independent datasets."""
import csv
from decimal import Decimal, InvalidOperation
import hashlib


def feature_id(label):
    """Remove vsearch size annotations and FASTA descriptions from local labels."""
    tokens = label.split()
    if not tokens or not tokens[0].split(';')[0]:
        raise ValueError('Feature labels must not be empty')
    return tokens[0].split(';')[0]


def sequence_ids(fasta):
    labels, sequences = {}, {}
    name, parts = None, []

    def store():
        if name is None:
            return
        sequence = ''.join(parts).upper()
        if not sequence:
            raise ValueError(f'Empty sequence for FASTA feature {name}')
        stable = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if name in labels and labels[name] != stable:
            raise ValueError(f'FASTA feature {name} identifies different sequences')
        labels[name] = stable
        sequences[stable] = sequence

    with open(fasta) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                store()
                name, parts = feature_id(line[1:]), []
            else:
                if name is None:
                    raise ValueError('FASTA sequence appears before its feature header')
                parts.append(line)
    store()
    if not labels:
        raise ValueError(f'No feature sequences found in {fasta}')
    return labels, sequences


def rewrite_features(fasta, table, output_fasta, output_table):
    """Rewrite a FASTA/table together, summing only identical nucleotide sequences."""
    labels, sequences = sequence_ids(fasta)
    samples, counts = None, {}
    with open(table, newline='') as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row:
                continue
            if samples is None:
                if row[0].startswith('#') and row[0] not in ('#OTU ID', '#OTU'):
                    continue
                samples = row[1:]
                if not samples or len(set(samples)) != len(samples):
                    raise ValueError('Feature table requires distinct sample columns')
                continue
            if row[0].startswith('#'):
                continue
            local_id = feature_id(row[0])
            if local_id not in labels:
                raise ValueError(f'Feature table ID {local_id} is missing from representative FASTA')
            if len(row) != len(samples) + 1:
                raise ValueError(f'Feature table row {local_id} has the wrong number of sample counts')
            try:
                values = [Decimal(value) for value in row[1:]]
            except InvalidOperation as error:
                raise ValueError(f'Invalid count in feature table row {local_id}') from error
            if any(not value.is_finite() or value < 0 or value != value.to_integral_value() for value in values):
                raise ValueError(f'Feature table row {local_id} requires non-negative integer counts')
            stable = labels[local_id]
            totals = counts.setdefault(stable, [0] * len(samples))
            for index, value in enumerate(values):
                totals[index] += int(value)
    if not counts:
        raise ValueError('Feature table contains no feature rows')
    with open(output_fasta, 'w') as handle:
        for stable in sorted(counts):
            handle.write(f'>{stable}\n{sequences[stable]}\n')
    with open(output_table, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['#OTU ID', *samples])
        for stable in sorted(counts):
            writer.writerow([stable, *counts[stable]])
    return {'feature_ids': set(counts), 'sample_ids': set(samples),
            'total_reads': sum(sum(values) for values in counts.values())}
