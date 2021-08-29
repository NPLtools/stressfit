echo off
rem	make package and install locally
rem	from Anaconda console, stressfit root, execute

rem	source package
python setup.py sdist
rem	binary package
python setup.py bdist_egg
rem install editable package
python -m pip install -e .
