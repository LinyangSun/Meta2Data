"""Offline XML parsing with Biopython and an unwritable container HOME."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIO_AVAILABLE = importlib.util.find_spec('Bio') is not None
APP_AVAILABLE = all(importlib.util.find_spec(name) is not None
                    for name in ('Bio', 'pandas', 'numpy', 'requests'))

# Use Biopython's bundled esearch.dtd. Any attempt to fetch a DTD from the
# internet in the child is rejected, so these regressions remain fully offline.
XML = b'''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN"
 "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult><Count>1</Count><RetMax>1</RetMax><RetStart>0</RetStart>
<IdList><Id>12345</Id></IdList><TranslationSet/>
<QueryTranslation>test</QueryTranslation></eSearchResult>
'''

CHILD = r'''
import io
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch
sys.path.insert(0, sys.argv[1])
from Bio import Entrez
assert 'Bio.Entrez.Parser' not in sys.modules, 'Parser already imported in fresh child'
before = dict(os.environ)
# /sys is a read-only kernel filesystem even for root in the test container.
try:
    Path(os.environ['HOME'], '.meta2data-cache-write-probe').mkdir()
except OSError:
    pass
else:
    Path(os.environ['HOME'], '.meta2data-cache-write-probe').rmdir()
    raise AssertionError('Regression requires an actually unwritable HOME')
entry = sys.argv[2]
if entry == 'pip':
    import py_16s
    initialize = py_16s._configure_entrez
elif entry == 'metadata':
    import metadata_downloader
    initialize = lambda: metadata_downloader.setup_entrez('offline-test-key')
else:
    from entrez_cache import configure_entrez_cache
    initialize = configure_entrez_cache
initialize()
from Bio.Entrez import Parser
first = Entrez.local_cache
assert first == Parser.DataHandler.directory
assert Path(first, 'Bio', 'Entrez', 'DTDs').is_dir()
assert Path(first, 'Bio', 'Entrez', 'XSDs').is_dir()
assert (Path(Entrez.__file__).parent / 'DTDs' / 'esearch.dtd').is_file()
with patch.object(Parser, 'urlopen', side_effect=AssertionError('Network access forbidden')):
    result = Entrez.read(io.BytesIO(bytes.fromhex(sys.argv[3])))
assert list(result['IdList']) == ['12345']
initialize()
assert Entrez.local_cache == first == Parser.DataHandler.directory
assert dict(os.environ) == before, 'Cache setup mutated global environment'
if entry == 'metadata':
    assert Entrez.api_key == 'offline-test-key'
print(json.dumps({'cache': first, 'home': os.environ['HOME']}))
'''


@unittest.skipUnless(BIO_AVAILABLE and Path('/sys').is_dir(),
                     'Requires Biopython and Linux read-only /sys')
class EntrezCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='entrez test ')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def parse_in_child(self, xdg=None, entry='helper'):
        env = os.environ.copy()
        env['HOME'] = '/sys'
        env['TMPDIR'] = str(self.root)
        env.pop('XDG_CACHE_HOME', None)
        if xdg is not None:
            env['XDG_CACHE_HOME'] = xdg
        result = subprocess.run(
            [sys.executable, '-c', CHILD, str(ROOT / 'scripts'), entry, XML.hex()],
            env=env, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads(result.stdout.strip().splitlines()[-1])

    def test_fresh_parser_honors_explicit_xdg_with_read_only_home(self):
        xdg = self.root / 'cache with spaces'
        result = self.parse_in_child(str(xdg))
        self.assertEqual(Path(result['cache']), xdg / 'meta2data' / 'biopython')
        self.assertEqual(result['home'], '/sys')

    def test_missing_xdg_uses_private_writable_fallback(self):
        result = self.parse_in_child()
        cache = Path(result['cache'])
        self.assertTrue(cache.is_relative_to(self.root))
        self.assertTrue(cache.parent.name.startswith('meta2data-entrez-'))
        # Temporary fallback directories are cleaned at interpreter exit.
        self.assertFalse(cache.exists())

    def test_read_only_xdg_uses_fallback(self):
        result = self.parse_in_child('/sys')
        self.assertTrue(Path(result['cache']).is_relative_to(self.root))
        self.assertFalse(Path(result['cache']).exists())

    def test_relative_xdg_is_not_interpreted_against_working_directory(self):
        result = self.parse_in_child('relative-cache-not-valid-xdg')
        self.assertTrue(Path(result['cache']).is_relative_to(self.root))

    def test_non_directory_xdg_uses_fallback(self):
        blocked = self.root / 'cache-file'
        blocked.write_text('keep this file')
        result = self.parse_in_child(str(blocked))
        self.assertTrue(Path(result['cache']).is_relative_to(self.root))
        self.assertEqual(blocked.read_text(), 'keep this file')

    @unittest.skipUnless(APP_AVAILABLE, 'Requires application dependencies')
    def test_public_platform_entrypoint_configures_cache_before_parser_import(self):
        xdg = self.root / 'pip-cache'
        result = self.parse_in_child(str(xdg), 'pip')
        self.assertEqual(Path(result['cache']), xdg / 'meta2data' / 'biopython')

    @unittest.skipUnless(APP_AVAILABLE, 'Requires application dependencies')
    def test_metadata_entrypoint_configures_cache_without_losing_api_key(self):
        xdg = self.root / 'metadata-cache'
        result = self.parse_in_child(str(xdg), 'metadata')
        self.assertEqual(Path(result['cache']), xdg / 'meta2data' / 'biopython')

    def test_loaded_parser_is_reconfigured_when_explicit_cache_changes(self):
        child = r'''
import io
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch
sys.path.insert(0, sys.argv[1])
from entrez_cache import configure_entrez_cache
from Bio import Entrez
first = configure_entrez_cache()
from Bio.Entrez import Parser
assert Parser.DataHandler.directory == first
os.environ['XDG_CACHE_HOME'] = sys.argv[2]
before = dict(os.environ)
second = configure_entrez_cache()
assert second != first
assert Parser.DataHandler.directory == Entrez.local_cache == second
assert dict(os.environ) == before
with patch.object(Parser, 'urlopen', side_effect=AssertionError('Network access forbidden')):
    assert list(Entrez.read(io.BytesIO(bytes.fromhex(sys.argv[3])))['IdList']) == ['12345']
print(json.dumps({'cache': second}))
'''
        env = {**os.environ, 'HOME': '/sys',
               'XDG_CACHE_HOME': str(self.root / 'first'), 'TMPDIR': str(self.root)}
        second = self.root / 'second'
        result = subprocess.run(
            [sys.executable, '-c', child, str(ROOT / 'scripts'), str(second), XML.hex()],
            env=env, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(Path(json.loads(result.stdout)['cache']),
                         second / 'meta2data' / 'biopython')


if __name__ == '__main__':
    unittest.main()
