# -*- coding: utf-8 -*-
"""
Shape definition for a sphere
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
from shapes import ShapeAbstract
import numpy as np
from numpy.linalg import norm
import shapes.tracing as tr


class ShapeSph(ShapeAbstract):
    shape_type = 'sphere'

    def __init__(self, radius):
        super().__init__()
        self.rad = radius

    def depthLocal(self, r):
        rn = norm(r, axis=1)
        res = self.rad - rn
        ins = np.array(res > 0, dtype=int)
        return [res, rn, ins]

    def cross(self, r, k):
        [tin, tout] = tr.crossSphere(self.rad, r, k)
        n = r.shape[0]
        tmp = np.array([tin, tout]).reshape((n, 2))
        return list(tmp)

    def rayPaths(self, r, ki, kf):
        [ti, ti2] = tr.crossSphere(self.rad, r, ki)
        [tf1, tf] = tr.crossSphere(self.rad, r, kf)
        sg1 = np.array(ti < 0, dtype=int)
        sg2 = np.array(tf > 0, dtype=int)
        path1 = -sg1*ti*norm(ki)
        path2 = sg2*tf*norm(kf)
        return [path1, path2]

    def plotContours(self, plt, proj, color, linestyle):
        theta = np.arange(0, 2*np.pi, 0.01)
        x = self.rad*np.cos(theta)
        y = self.rad*np.sin(theta)
        plt.plot(x, y, color=color, linestyle=linestyle)

