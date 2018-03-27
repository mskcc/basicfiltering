# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='basicfiltering',
    version='0.2.0',
    description='Basic false-positive filters for VCFs/TXTs from somatic variant callers',
    long_description=readme,
    author='Ronak Shah, Cyriac Kandoth',
    author_email='ckandoth@gmail.com',
    url='https://github.com/mskcc/basicfiltering',
    license=license,
    install_requires=['nose==1.3.7', 'pyvcf==0.6.8', 'pandas==0.19.2'],
    packages=find_packages(exclude=('tests', 'docs'))
)
