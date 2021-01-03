# -*- coding: utf-8 -*-

import stressfit.matrix.components as c
import numpy as np
import pytest


def genRay(x=1.5, z=-2.0, a=0.001, dll=0.01, time=10.0, v0=2.0):
    """Generate test trajectory at the sample."""    
    X = [x, z, a, dll, time]
    return np.array(X)


def ans_slit(X, L=1000, v0=2.0, alpha=60*c.deg, omega=4*c.deg):
    alpha2 = 0.5*alpha - omega
    ca = np.cos(alpha2)
    sa = np.sin(alpha2)
    tt = np.tan(0.5*alpha)
    x0 = sa*X[0] + ca*X[1]
    z0 = -ca*X[0] + sa*X[1] 
    x1 = x0 + X[2]*L + 2*L*tt*X[3]
    a1 = X[2] + 2*tt*X[3]
    t1 = X[4] + L/v0*X[3] - z0/v0
    ans = [x1, a1, X[3], t1]
    return np.array(ans)
    

def test_slit(L=1000, v0=2.0, alpha=60*c.deg, omega=4*c.deg):
    comp = c.Slit(distance=L, width=2)
    alpha2 = 0.5*alpha - omega
    C = c.CD(0, alpha2, 0.5*alpha, v0)
    comp.initialize(C, v0)
    X = genRay(v0=v0)
    X1 = comp.Cout.dot(X)
    ans = ans_slit(X, L=L, v0=v0, alpha=alpha, omega=omega)
    assert X1 == pytest.approx(ans)
    







