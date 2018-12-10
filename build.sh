#!/bin/bash

rm -r fastatools.egg-info/
python setup.py clean --all
python setup.py sdist --verbose

