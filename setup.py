from setuptools import setup, find_packages

setup(
    name='ot_service_wrapper',
    version='0.1',
    packages=find_packages(exclude=['tests*']),
    license='gpl3',
    description='Python wrapper for Open Tree services',
    author='CJ Grady',
    author_email='cjgrady@ku.edu'
)
