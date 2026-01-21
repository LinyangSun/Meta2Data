#!/usr/bin/env python
"""Setup script for Meta2Data package.

This setup.py is needed to install bash scripts from the bin/ directory,
as pyproject.toml doesn't support non-Python executables well.
"""
from setuptools import setup

setup(
    scripts=[
        'bin/Meta2Data',
        'bin/Meta2Data-AmpliconPIP',
        'bin/Meta2Data-MetaDL',
    ],
)
