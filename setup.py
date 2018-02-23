#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
from pip.req import parse_requirements

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

install_requirements = parse_requirements('requirements.txt')
test_requirements = parse_requirements('requirements-dev.txt')

setup(
    author="Ingo Scholtes",
    author_email='ischoltes@ethz.ch',
    license='AGPL-3.0+',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: AGPL-3.0+ License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="A python package for the analysis of sequential data on pathways and "
                "temporal networks from the perspective of higher-order network models.",
    install_requires=install_requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    python_requires='>=3.5',
    keywords='network analysis temporal networks pathways sequence modeling graph mining',
    name='pathpy',
    packages=find_packages(include=['pathpy']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/IngoScholtes/pathpy',
    version='1.1.2',
    zip_safe=False,
)
