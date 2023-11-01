# -*- coding: utf-8 -*-
"""
Module for optional binding with CUDA GPU.

Handles switching between numpy and cupy modes, with fallback to numpy if 
CUDA GPU and cupy is not present. 

Created on Mon Sep 18 12:59:48 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""
print('importing cuda')
from . import options

if options['gpu']:
    try:
        import cupy as cp
        _gpu = True
        print('Imported CuPy.')
    except:
        import numpy as cp
        options['gpu'] = False
        print('Cannot import CuPy, switched to NumPy.') 
else:
    import numpy as cp
    options['gpu'] = False
    print('Imported NumPy.')

def has_gpu():
    """Return True if cuda and cupy is available."""
    return options['gpu']

def asnumpy(a):
    """Convert a cupy array to numpy."""
    if options['gpu']:
        return cp.asnumpy(a)
    else:
        return a

def to_numpy(a):
    """Convert a list or dict of cupy arrays to numpy."""
    if options['gpu']:
        if isinstance(a,list):
            res = []
            for x in a:
                res.append(asnumpy(x))
        elif isinstance(a,dict):
            res = {}
            for key in a:
                res[key] = asnumpy(a[key])
        return res
    else:
        return a
