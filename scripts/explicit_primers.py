"""Trim user-supplied primers using the same file layout as automatic detection."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from read_layout import discover, layout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--forward', required=True)
    p.add_argument('--reverse', default='')
    p.add_argument('--detect-only', action='store_true')
    p.add_argument('--mixed-orientation', action='store_true')
    a = p.parse_args()
    rows = discover(a.input)
    mode = layout(rows)
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    reverse_used = bool(a.reverse) and (mode == 'PE' or a.mixed_orientation)

    def entry(sequence, used):
        applied = used and not a.detect_only
        return dict(name='explicit' if sequence else 'none', consensus=sequence,
                    length=len(sequence), detected=bool(sequence),
                    trim_length=None if applied else 0,
                    trim_method='cutadapt_variable' if applied else 'none')

    info = dict(mode=mode, layout=mode, status='completed', source='explicit',
                reason='Explicit primers recorded; reads copied unchanged' if a.detect_only
                else 'Explicit primers applied with cutadapt; trim lengths vary by read',
                trim_method='none' if a.detect_only else 'cutadapt_variable',
                forward_primer=entry(a.forward, True),
                reverse_primer=entry(a.reverse, reverse_used),
                detect_only=a.detect_only)

    for row in rows:
        r1, r2 = row['r1'], row['r2']
        if a.detect_only:
            for src in (r1, r2):
                if src:
                    shutil.copyfile(src, out / Path(src).name)
            continue
        cmd = ['cutadapt', '-g', a.forward, '-o', str(out / Path(r1).name)]
        if r2:
            if a.reverse:
                cmd += ['-G', a.reverse]
            cmd += ['-p', str(out / Path(r2).name), r1, r2]
        else:
            if a.mixed_orientation and a.reverse:
                cmd += ['-g', a.reverse]
            cmd += [r1]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as error:
            reason = f'cutadapt failed for {Path(r1).name} (exit code {error.returncode})'
            info.update(status='failed', reason=reason)
            (out / 'primer_info.json').write_text(json.dumps(info, indent=2) + '\n')
            print(f'Error: {reason}', file=sys.stderr)
            return error.returncode
    (out / 'primer_info.json').write_text(json.dumps(info, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
