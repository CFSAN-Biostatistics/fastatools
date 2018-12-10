#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'Biopython',
]

test_requirements = [
    "pytest",
]

setup(
    name='fastatools',
    version='1.2.0',
    description="Tools for working with fasta files.",
    long_description=readme + '\n\n' + history,
    author="Steve Davis",
    author_email='steven.davis@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/fastatools',
    packages=[
        'fastatools',
    ],
    package_dir={'fastatools':
                 'fastatools'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['bioinformatics', 'NGS', 'fastatools'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    entry_points={'console_scripts': ['fastatools = fastatools.cli:main']},
    setup_requires=["pytest-runner"],
    tests_require=test_requirements
)
