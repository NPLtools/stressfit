# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:23:11 2020
@author: Jan Saroun, saroun@ujf.cas.cz
"""

from setuptools import setup, find_packages

def get_version():
    """Read version from package resources."""
    try:
        import json
        fn = r'stressfit/resources/conf/setup.json'
        with open(fn, 'r') as f:
            lines = f.readlines() 
        res = json.loads('\n'.join(lines))
        version = res['version']
    except:
        version = 'unknown'
    return version

stressfit_data = ['resources/*/*']

other_data = ['stressfit_example1.py', 
              'LICENSE', 
              'README.md', 
              'CHANGELOG']

requires = [
'numpy', 
'matplotlib',
'lmfit',
'pathlib',
'colorama',
'ipywidgets',
'ipysheet'
]

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

version_string = get_version()

setup(
    name = 'stressfit',
    version = version_string,
    description = 'Fitting of residual stress distributions measured by neutron diffraction.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    author = 'Jan Saroun',
    author_email = 'saroun@ujf.cas.cz',
    license = 'MPL-2.0',
    url = 'https://github.com/NPLtools/stressfit',
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
    ],
    platforms = ['any'],
    keywords = 'residual stress, neutron, pseudo-strain, surface effect', 
    python_requires = '>=3.8',
    install_requires = requires,
    packages = find_packages(),
    package_data = {'stressfit': stressfit_data}, 
    data_files = [('.',other_data)]
)
