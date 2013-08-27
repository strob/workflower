#!/usr/bin/env python

import distutils.core

distutils.core.setup(
    name='workflower',
    version='0.1',
    author='Robert M Ochshorn',
    author_email='rmo@numm.org',
    packages=['workflower'],
    package_data = {"workflower": ["data/index.html"]},
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        ],
    license='GPL')
