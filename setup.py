#!/usr/bin/env python

import distutils.core
import setuptools
import glob

def get_datalist():
    return [ X.replace("workflower/", "") # relative paths
             for X in 
             glob.glob("workflower/data/*") ]

distutils.core.setup(
    name='workflower',
    version='0.1',
    author='Robert M Ochshorn',
    author_email='rmo@numm.org',
    packages=['workflower'],
    package_data = {"workflower": get_datalist()},
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        ],
    license='GPL')
