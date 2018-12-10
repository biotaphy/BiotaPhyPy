# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ot_service_wrapper',
    version='0.1',
    description='Biotaphy package for wrapping Open Tree services',
    long_description=readme,
    author='Biotaphy Team',
    author_email='cjgrady@ku.edu',
    url='https://github.com/biotaphy/ot_service_wrapper',
    license=license,
    packages=find_packages(exclude=('tests')),
    install_requires=[
        'dendropy>=4.0.0'
    ]
)

