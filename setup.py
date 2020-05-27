# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

with open('LICENSE', 'r', encoding='utf-8') as f:
    license = f.read()

setup(
    name='biotaphy',
    version='1.2.0',
    description='BiotaPhy Python Package',
    long_description=readme,
    author='Biotaphy Team',
    author_email='cjgrady@ku.edu',
    url='https://github.com/biotaphy/BiotaPhyPy',
    license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample_data')),
    scripts=['bin/ancestral_distribution.py', 'bin/phylo_beta_diversity.py'],
    install_requires=[
        'dendropy>=4.0.0',
        'matplotlib',
        'numpy>=1.11.0',
        'scipy>=1.0.0']
)
