# -*- coding: utf-8 -*-
"""Setup module for BiotaPhyPy."""
from setuptools import setup, find_packages
import versioneer


with open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

with open('LICENSE', 'r', encoding='utf-8') as f:
    license = f.read()

setup(
    name='biotaphy',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='BiotaPhy Python Package',
    long_description=readme,
    author='Biotaphy Team',
    author_email='cjgrady@ku.edu',
    url='https://github.com/biotaphy/BiotaPhyPy',
    license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample_data')),
    install_requires=[
        'defusedxml',
        'gdal',
        'dendropy>=4.0.0',
        'matplotlib',
        'numpy>=1.11.0',
        'scipy>=1.0.0',
        'rtree',
        'specify-lmpy>=3.1.1',
    ]
)
