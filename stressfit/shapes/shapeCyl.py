# -*- coding: utf-8 -*-
"""
Shape definition for a cylinder (y-axis vertical)
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapeCyl(ShapeAbstract):
    """
    Shape definition for a cylinder (y-axis vertical).
    """
    
    shape_type = 'cylinder'

    def __init__(self, radius, height):
        super().__init__()
        self.rad = radius
        self.height = height

    def depthLocal(self, r):
        r1 = r[:, 0::2]
        rn = norm(r1, axis=1)
        res = self.rad - rn
        ins = np.array(res > 0, dtype=int)
        return [res, rn, ins]

    def cross(self, r, k):
        [tin, tout] = tr.crossCylinder(self.rad, r, k)
        # apply other boundaries
        ty = tr.crossLayer(self.height, r, k, 1)
        tin = np.maximum(tin, ty[0])
        tout = np.minimum(tout, ty[1])
        n = r.shape[0]
        tmp = np.array([tin, tout]).reshape((n, 2))
        return list(tmp)

    def rayPaths(self, r, ki, kf):
        [ti1, ti2] = tr.crossCylinder(self.rad, r, ki)
        tiy = tr.crossLayer(self.height, r, ki, 1)
        tin = np.maximum(ti1, tiy[0].T).T
        [tf1, tf2] = tr.crossCylinder(self.rad, r, kf)
        tfy = tr.crossLayer(self.height, r, kf, 1)
        tout = np.minimum(tf2, tfy[1].T).T
        sg1 = np.array(tin < 0, dtype=int)
        sg2 = np.array(tout > 0, dtype=int)
        path1 = -sg1*tin*norm(ki)
        path2 = sg2*tout*norm(kf)
        return [path1, path2]

    def plotContours2(self, plt, proj, color, linestyle):
        theta = np.arange(0, 2*np.pi, 0.01)
        if ((proj == 0) or (proj == 2)):
            wx = self.rad
            hx = self.height*0.5
            x = [-wx, -wx, wx, wx, -wx]
            y = [-hx, hx, hx, -hx, -hx]
            plt.plot(x, y, color=color, linestyle=linestyle)
        elif (proj == 1):
            x = self.rad*np.cos(theta)
            y = self.rad*np.sin(theta)
            plt.plot(x, y, color=color, linestyle=linestyle)

    def plotContours(self, ax, proj, color, linestyle):
        theta = np.arange(-np.pi/2, np.pi/2, 0.01)
        gray = (0.5, 0.5, 0.5, 0.15)
        if ((proj == 0) or (proj == 2)):
            w = self.rad
            h = self.height*0.5
            h = self.height*0.5           
            # fill gray
            x = [-w, -w, w, w]
            y = [-h, h, h, -h]
            ax.fill_between([-w, w], [-h, -h], [h, h], facecolor=gray) 
            # contour
            x = [-w, -w, w, w, -w]
            y = [-h, h, h, -h, -h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='a')          
        elif (proj == 1):
            x = self.rad*np.sin(theta)
            y2 = self.rad*np.cos(theta)
            y1 = -y2
            ax.fill_between(x, y1, y2, facecolor=gray)
            ax.plot(x, y1, color=color, linestyle=linestyle, label='a')
            ax.plot(x, y2, color=color, linestyle=linestyle, label='b')
            