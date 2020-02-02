STRESSFIT library
Written by: Jan Saroun, Nuclear Physics Institute CAS, Rez, saroun@ujf.cas.cz
Date: 11/09/2017 
Last update: 22/5/2018

==========
Contents
==========

Execution scripts:
------------------
stressfit-test.py:	MC convolution demonstrator
stressfit-worksheet.py:  Template script for data fitting
available also as Jupyter notebook  files, *.ipynb (see http://jupyter.org)

Libraries:
-----------
sample.py: handling of sample properties and MC convolution
mccfit.py: classes and functions for data fitting using MC convolution
graphs.py: some plot functions and definitions
./shapes/*.py: classes with sample shape definitions and methods

Example data:
-------------
./input/events_EX_1mm.dat: sampling events for ENGIN-X, Fe(211), gauge volume ~ 1x1x5 mm 
./input/events_SS_1mm.dat: sampling events for STRESS-SPEC, Fe(211), lambda=1.68 A, gauge volume  ~ 1x1x5 mm 
Fe_mu.dat: Total removal cross-section table for alpha-Fe
./input/*.dat: example input data (synthetic data simulated by SIMRES)

Other:
-------
fit_centre.py:	example script for fitting of sample centre from two perpendicular intensity scans.
eps2sig.py: script for conversion of strain distributions to stresses
./input/strain_table.dat: strain and intensity distribution table used to simulate synthetic data in SIMRES
