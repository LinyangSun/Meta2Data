#!/usr/bin/env python
"""Setup script for Meta2Data package.

This setup.py is needed to install bash scripts from the bin/ directory
and data files from the test/ and docs/ directories, as pyproject.toml
doesn't support non-Python executables and data files well.
"""
from setuptools import setup
import os
import glob

# Get the directory containing this setup.py
here = os.path.abspath(os.path.dirname(__file__))

# Prepare data files
data_files = []

# Test CSV file
test_csv = os.path.join(here, 'test', 'ampliconpiptest.csv')
if os.path.exists(test_csv):
    data_files.append(('share/Meta2Data/test', ['test/ampliconpiptest.csv']))

# Docs directory (includes reference sequences and documentation)
docs_files = glob.glob('docs/*')
if docs_files:
    data_files.append(('share/Meta2Data/docs', docs_files))

# Scripts directory (includes AmpliconFunction.sh, py_16s.py, run.sh, taxonomy.sh)
scripts_files = glob.glob('scripts/*.sh') + glob.glob('scripts/*.py')
if scripts_files:
    data_files.append(('share/Meta2Data/scripts', scripts_files))

setup(
    scripts=[
        'bin/Meta2Data',
        'bin/Meta2Data-AmpliconPIP',
        'bin/Meta2Data-ggCOMBO',
        'bin/Meta2Data-MetaDL',
    ],
    data_files=data_files,
)
