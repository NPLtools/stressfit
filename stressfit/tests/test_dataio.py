# -*- coding: utf-8 -*-
"""
Test module dataio.py.

"""

import stressfit.dataio as io
import numpy as np
import json
import pytest
from pathlib import Path
    

def test_load_data():
    """Read data using numpy with defined range of rows ansd columns."""
    data = io.load_data("int_B_axi.dat", kind='data', rows=[1,4],
                        usecols=(1,2))  
    assert data.shape == (4,2)


def test_data_rw(tmp_path):
    """Read and write of data using numpy."""
    d:Path = tmp_path 
    print('Temporary output path: '+d.as_posix())
    if not d.exists():
        d.mkdir()
    data = io.load_data("int_B_axi.dat", kind='data')
    io.save_data(data, 'tmp.dat', source=__file__, path=d)
    data1 = io.load_data("tmp.dat", path=d)
    # compare values
    diff = np.isclose(data, data1, rtol=1e-6, atol=1e-12)
    assert diff.all()


def test_table_rw(tmp_path):
    """Read and write of Table."""
    d:Path = tmp_path 
    print('Temporary output path: '+d.as_posix())
    if not d.exists():
        d.mkdir()
    data = io.load_table("test.tab", kind='tables')
    io.save_table(data, 'tmp.dat', source=__file__, path=d)
    data1 = io.load_table("tmp.dat", path=d)
    # check table size
    assert len(data.rows)>3
    assert data.size == data1.size
    # compare values
    diff = np.isclose(data.cells, data1.cells, rtol=1e-6, atol=1e-12)
    assert diff.all()
    # compare headers
    if len(data.colHeaders)>0:
        for i in range(len(data.colHeaders)):
            assert data.colHeaders[i] == data1.colHeaders[i]
    if len(data.rowHeaders)>0:
        for i in range(len(data.rowHeaders)):
            assert data.rowHeaders[i] == data1.rowHeaders[i]

def test_params_rw(tmp_path):
    """Read and write parameters."""
    d:Path = tmp_path 
    print('Temporary output path: '+d.as_posix())
    if not d.exists():
        d.mkdir()
    data = io.load_params("ENGINX.par", kind='instruments')
    io.save_params(data, 'tmp.dat', source=__file__, path=d)
    data1 = io.load_params("tmp.dat", path=d)
    for k in data:
        assert k in data1
        assert data[k].value == pytest.approx(data1[k].value)
        assert data[k].descr == data1[k].descr
    

def test_load_text():
    """Test JSON input by loading as a text and decoding."""
    path = io.__instruments
    txt = io.load_text("ENGINX.json", path=path)
    cont = "\n".join(txt)
    data = json.loads(cont)
    comp =  ['src', 'guide', 'col_0', 'col_1', 'col_2', 'det']
    for c in comp:
        assert c in data

# test_load_data()
# test_data_rw(Path.home())
# test_table_rw(Path.home())
# test_params_rw(Path.home())

