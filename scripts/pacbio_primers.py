"""Select the known/explicit primer required to orient DADA2 CCS reads."""
import json
from pathlib import Path
import sys
from entropy_primer_detect import PRIMERS_F, PRIMERS_R, reverse_complement_iupac


def sequence(entry, database):
    name = entry.get('name', '')
    if name == 'explicit':
        return entry['consensus']
    rc = name.endswith('_RC')
    original = name[:-3] if rc else name
    seq = dict(database).get(original, '')
    return reverse_complement_iupac(seq) if rc else seq


if __name__ == '__main__':
    info = json.loads(Path(sys.argv[1]).read_text())
    forward = sequence(info.get('forward_primer', {}), PRIMERS_F)
    reverse = sequence(info.get('reverse_primer', {}), PRIMERS_R)
    if not forward:
        sys.exit('SKIP: DADA2 CCS requires a known or explicit forward primer; supply --primer-fwd or use --vsearch.')
    # No reverse adapter is required by denoise-ccs when it was not detected.
    print(forward + '\t' + reverse)
