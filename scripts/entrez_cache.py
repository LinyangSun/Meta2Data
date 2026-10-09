"""Configure Biopython's XML cache without requiring a writable home directory.

Bio.Entrez.Parser reads Entrez.local_cache when its DataHandler class is first
imported; changing XDG variables alone does not redirect that cache. Keep all
changes local to Biopython and leave HOME and the process environment untouched.
"""

import os
from pathlib import Path
import tempfile
import threading


_lock = threading.Lock()
_fallback = None


def _prepare(directory):
    """Create and actually probe both directories used by Biopython's parser."""
    directory = Path(directory)
    for kind in ('DTDs', 'XSDs'):
        child = directory / 'Bio' / 'Entrez' / kind
        child.mkdir(parents=True, exist_ok=True)
        # os.access() alone cannot reliably identify read-only container binds.
        with tempfile.TemporaryFile(dir=child) as probe:
            probe.write(b'1')
            probe.flush()
    return str(directory.resolve())


def configure_entrez_cache():
    """Return a usable cache path and configure new or already-loaded parsers.

    Prefer an absolute XDG_CACHE_HOME. If it is absent, relative, or unwritable,
    use a private temporary directory retained for this process and cleaned on
    exit. The same fallback is reused on subsequent calls.
    """
    global _fallback
    with _lock:
        directory = None
        xdg = os.environ.get('XDG_CACHE_HOME')
        if xdg and Path(xdg).is_absolute():
            try:
                directory = _prepare(Path(xdg) / 'meta2data' / 'biopython')
            except OSError:
                # A user may bind an existing shared cache read-only. A private
                # temporary cache still permits bundled/offline DTD parsing.
                pass
        if directory is None:
            if _fallback is None:
                _fallback = tempfile.TemporaryDirectory(prefix='meta2data-entrez-')
            directory = _prepare(Path(_fallback.name) / 'biopython')

        from Bio import Entrez
        # This assignment MUST precede Parser import: its metaclass initializes
        # the directory and can raise EROFS before Entrez.read() parses any XML.
        Entrez.local_cache = directory
        from Bio.Entrez import Parser
        # Also update an already-imported parser; setting local_cache alone is
        # insufficient after DataHandler's metaclass has initialized.
        Parser.DataHandler.directory = directory
        return directory
