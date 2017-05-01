# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='basicfiltering',
    version='0.1.5',
    description='Package for filtering multiple variant calling tools',
    long_description=readme,
    author='Ronak Shah',
    author_email='rons.shah@gmail.com',
    url='https://github.com/rhshah/basicfiltering',
    license=license,
    install_requires=['nose==1.3.7', 'pyvcf==0.6.7', 'pandas==0.16.2', 'coloredlogs==5.2','codecov==2.0.5', 'coverage==4.3.4'],
    packages=find_packages(exclude=('tests', 'docs'))
)
