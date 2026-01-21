#!/usr/bin/env python
"""Setup script for Meta2Data package.

This setup.py is needed to install bash scripts from the bin/ directory
and data files from the test/ directory, as pyproject.toml doesn't support
non-Python executables and data files well.
"""
from setuptools import setup
import os

# Get the directory containing this setup.py
here = os.path.abspath(os.path.dirname(__file__))

# Read test CSV file path
test_csv = os.path.join(here, 'test', 'ampliconpiptest.csv')
test_data = []
if os.path.exists(test_csv):
    test_data = [('share/Meta2Data/test', ['test/ampliconpiptest.csv'])]

setup(
    scripts=[
        'bin/Meta2Data',
        'bin/Meta2Data-AmpliconPIP',
        'bin/Meta2Data-MetaDL',
    ],
    data_files=test_data,
)
