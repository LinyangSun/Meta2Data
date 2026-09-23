"""Download only the GG2 16S backbone needed before guarded AmpliconPIP."""
import argparse
import fcntl
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile


def valid(path):
    try:
        with zipfile.ZipFile(path) as archive:
            metadata = [n for n in archive.namelist() if n.count('/') == 1 and n.endswith('/metadata.yaml')]
            return (len(metadata) == 1 and b'FeatureData[Sequence]' in archive.read(metadata[0])
                    and archive.testzip() is None)
    except (OSError, zipfile.BadZipFile):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', required=True)
    args = parser.parse_args()
    directory = Path(args.db).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / '2024.09.backbone.full-length.fna.qza'
    with (directory / '.adapter_reference.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if valid(target):
            return
        with tempfile.NamedTemporaryFile(dir=directory, suffix='.qza', delete=False) as stream:
            temporary = Path(stream.name)
        try:
            url = 'https://ftp.microbio.me/greengenes_release/2024.09/' + target.name
            subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--output', str(temporary), url], check=True)
            if not valid(temporary):
                raise ValueError('Downloaded GG2 backbone is incomplete or has the wrong artifact type')
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
