# -*- coding: utf-8 -*-
"""
Test module dataio.py.

"""

import stressfit.dataio as io
import numpy as np
import pytest


def test_load_data():
    data = io.load_data("int_B_axi.dat", kind='input', rmin=1, rmax=4, 
                       usecols=(1,2))  
    assert data.shape == (4,2)

def test_read_dict():
    data = io.read_dict("ENGINX.par", kind='instrument')
    assert len(data.keys()) == 19

def test_read_table():
    data = io.read_table("Fe_mu.dat", kind='table')
    assert len(data.colHeaders) == 2

