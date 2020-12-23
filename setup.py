# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:23:11 2020
@author: Jan Saroun, saroun@ujf.cas.cz
"""

from setuptools import setup, find_packages

stressfit_data = ['tables/*', 'examples/*']

other_data = ['stressfit_example1.py', 
                'LICENSE', 
                'README.md', 
                'CHANGELOG']

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='stressfit',
    version='1.0.0',
    description='Fitting of residual stress distributions measured by neutron diffraction.',
    long_description = long_description,
    long_description_content_type='text/markdown',
    author='Jan Saroun',
    author_email='saroun@ujf.cas.cz',
    license='MPL-2.0',
    url='https://github.com/NPLtools/stressfit',
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
    ],
    platforms=['any'],
    keywords='residual stress, neutron, pseudo-strain, surface effect', 
    python_requires='>=3.6',
    install_requires=['numpy', 'matplotlib'],
    packages=find_packages(),
    package_data={'stressfit': stressfit_data}, 
    data_files=[('.',other_data)]
)
