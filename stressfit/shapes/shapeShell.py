# -*- coding: utf-8 -*-
"""
Shape definition for a hollow sphere
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapeShell(ShapeAbstract):
    """Define concentric hollow sphere.
    
    Parameters
    ----------
    Rin: float
        Inner (hole) radius [mm].
        
    Rout: float
        Outer radius [mm].
    """
    shape_type = 'shell'

    def __init__(self, Rin=4.0, Rout=8.0):
        super().__init__()
        self.R1 = Rin
        self.R2 = Rout

    def update(self, **kwargs):
        """Update parameters."""
        if 'Rin' in kwargs:
            self.R1 = kwargs['Rin']
        if 'Rout' in kwargs:
            self.R2 = kwargs['Rout']
            
    def depthLocal(self, r):
        rn = norm(r, axis=1)
        r1 = self.R1 - rn
        r2 = self.R2 - rn
        ins = np.array((r2 > 0) & (r1 < 0), dtype=int)
        return [r2, r1, ins]

    def cross(self, r, k):
        n = r.shape[0]
        [tin, tout] = tr.crossShell(self.R1, self.R2, r, k)
        cm = []
        for i in range(n):
            tmp = [[tin[i, 0], tout[i, 0]], [tin[i, 1], tout[i, 1]]]
            cm.append(np.array(tmp))
        return cm

    def rayPaths(self, r, ki, kf):
        [ti1, ti2] = tr.crossShell(self.R1, self.R2, r, ki)
        [tf1, tf2] = tr.crossShell(self.R1, self.R2, r, kf)
        path1 = np.minimum(ti2, 0) - ti1
        path2 = tf2 - np.maximum(tf1, 0)
        mi = np.array(path1 > 0, dtype=int)
        mf = np.array(path2 > 0, dtype=int)
        res1 = np.sum(path1*mi, axis=1)*norm(ki)
        res2 = np.sum(path2*mf, axis=1)*norm(kf)
        return [res1, res2]

    def plotContours(self, plt, proj, color, linestyle):
        theta = np.arange(0, 2*np.pi, 0.01)
        x = self.R1*np.cos(theta)
        y = self.R1*np.sin(theta)
        plt.plot(x, y, color=color, linestyle=linestyle, label='a')
        x = self.R2*np.cos(theta)
        y = self.R2*np.sin(theta)
        plt.plot(x, y, color=color, linestyle=linestyle, label='b')

