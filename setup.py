#!/usr/bin/env python
from setuptools import setup
import glob

scripts = glob.glob("*.p*")

setup(
    name='RaGOO',
    version='v1.1',
    description='A tool to order and orient genome assembly contigs via minimap2 alignments to a reference genome.',
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    packages=['ragoo_utilities'],
    package_dir={'ragoo_utilities': 'ragoo_utilities/'},
    install_requires=[
              'intervaltree',
              'numpy',
          ],
    scripts=scripts,
    zip_safe=True
)
