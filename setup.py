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
    # TODO: put package test requirements here
]

setup(
    name='fastatools',
    version='0.1.0',
    description="Tools to manipulate fasta files.",
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
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    entry_points={'console_scripts': ['fastatools = fastatools.fastatools:main']},
    test_suite='tests',
    tests_require=test_requirements
)
