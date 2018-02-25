#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file, open('HISTORY.rst') as history_file:
    readme = readme_file.read()
    history = history_file.read()

install_requirements = ['numpy', 'scipy']

setup_requirements = ['pytest-runner']

setup(
    author="Ingo Scholtes",
    author_email='ischoltes@ethz.ch',
    license='AGPL-3.0+',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="A python package for the analysis of sequential data on pathways and "
                "temporal networks from the perspective of higher-order network models.",
    install_requires=install_requirements,
    setup_requires=setup_requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    python_requires='>=3.5',
    keywords='network analysis temporal networks pathways sequence modeling graph mining',
    name='pathpy',
    packages=find_packages(include=['pathpy']),
    test_suite='tests',
    url='https://github.com/IngoScholtes/pathpy',
    version='1.1.1',
    zip_safe=False
)
