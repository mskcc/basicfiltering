# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='basicfiltering',
    version='0.1.0',
    description='Package for filtering multiple variant calling tools',
    long_description=readme,
    author='Ronak Shah',
    author_email='rons.shah@gmail.com',
    url='https://github.com/rhshah/basicfiltering',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)