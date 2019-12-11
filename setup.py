# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='biotaphy',
    version='2.0.0',
    description='Biotaphy package for various computations',
    long_description=readme,
    author='Biotaphy Team',
    author_email='cjgrady@ku.edu',
    url='https://github.com/biotaphy/analyses',
    license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample_data')),
    scripts=['bin/ancestral_distribution.py'],
    install_requires=[
        'dendropy>=4.0.0',
        'matplotlib',
        'numpy>=1.11.0',
        'scipy>=1.0.0']
)
