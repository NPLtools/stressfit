# -*- coding: utf-8 -*-
"""
Shape definition for a planparallel plate
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
from shapes import ShapeAbstract
import numpy as np
from numpy.linalg import norm
import shapes.tracing as tr

_eps = 1.0e-12


class ShapePlate(ShapeAbstract):
    shape_type = 'plate'

    def __init__(self, thickness):
        super().__init__()
        self.t = thickness

    def depthLocal(self, r):
        d1 = 0.5*self.t-r[:, 2]
        d2 = -0.5*self.t-r[:, 2]
        ins = np.array((d1 > 0) & (d2 < 0), dtype=int)
        return [d1, d2, ins]

    def cross(self, r, k):
        n = r.shape[0]
        t = tr.crossLayer(self.t, r, k, 2)
        tmp = np.array(t).reshape((n, 2))
        return list(tmp)

    def rayPaths(self, r, ki, kf):
        [ti, ti2] = tr.crossLayer(self.t, r, ki, 2)
        [tf1, tf] = tr.crossLayer(self.t, r, kf, 2)
        sg1 = np.array(ti < 0, dtype=int)
        sg2 = np.array(tf > 0, dtype=int)
        path1 = -sg1*ti*norm(ki)
        path2 = sg2*tf*norm(kf)
        return [path1, path2]

    
    def plotContours2(self, plt, proj, color, linestyle):
        if (proj==2):
            return
        ax = plt.axes()
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        
        if (proj == 0):
            x = np.array([-0.5*self.t, -0.5*self.t])
            y = np.array([ylim[0], ylim[1]])
            plt.plot(x, y,color=color, linestyle=linestyle, marker=None, label='a')
            plt.plot(-x, y,color=color, linestyle=linestyle, marker=None, label='b')
        else:
            x = np.array([xlim[0], xlim[1]])
            y = np.array([-0.5*self.t, -0.5*self.t])
            plt.plot(x, y,color=color, linestyle=linestyle, marker=None, label='a')
            plt.plot(x, -y,color=color, linestyle=linestyle, marker=None, label='b')


    def plotContours(self, ax, proj, color, linestyle):
        if (proj==2):
            return
        gray = (0.2, 0.2, 0.2, 0.15)
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        if (proj == 0):
            x = np.array([-0.5*self.t, -0.5*self.t])
            y = np.array([ylim[0], ylim[1]])
            # fill gray
            ax.fill_between([-0.5*self.t,0.5*self.t], [ylim[0], ylim[0]], [ylim[1], ylim[1]], facecolor=gray) 
            ax.plot(x, y,color=color, linestyle=linestyle, marker=None, label='a')
            ax.plot(-x, y,color=color, linestyle=linestyle, marker=None, label='b')
        else:
            x = np.array([xlim[0], xlim[1]])
            y = np.array([-0.5*self.t, -0.5*self.t])
            # fill gray
            ax.fill_between([xlim[0], xlim[1]], [-0.5*self.t, -0.5*self.t], [0.5*self.t, 0.5*self.t], facecolor=gray) 
            ax.plot(x, y,color=color, linestyle=linestyle, marker=None, label='a')
            ax.plot(x, -y,color=color, linestyle=linestyle, marker=None, label='b')


        