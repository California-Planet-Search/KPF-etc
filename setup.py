# -*- coding: utf-8 -*-
import logging
import sys
from io import open
from os import path

try:
    from setuptools import setup, find_packages
except ImportError:
    logging.exception('Please install or upgrade setuptools or pip')
    sys.exit(1)

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='kpf_etc',
    version='1.0',
    description='A package for estimating exposure times for the Keck Planet Finder',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/California-Planet-Search/KPF-etc',
    author='Sam Halverson',
    author_email='samuel.halverson@jpl.nasa.gov',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.7',
    ],
    packages=['kpf_etc',
              'kpf_etc.grids'],
    
    package_data={
      # If any package contains *.txt, *.rst or *.fits files, include them:
      '': ['*.txt', '*.yaml', '*.fits', '*.pdf']
    },
    include_package_data=True,

    python_requires='>=2.7',
    install_requires=[
        'astropy',
        'matplotlib',
        'numpy',
        'scipy'
    ]
)
