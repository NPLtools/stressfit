# -*- coding: utf-8 -*-
"""
Shape definition for a hollow cylinder
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
from shapes import ShapeAbstract
import numpy as np
from numpy.linalg import norm
import shapes.tracing as tr


class ShapeShellCyl(ShapeAbstract):
    shape_type = 'shell'

    def __init__(self, Rin, Rout, height):
        super().__init__()
        self.R1 = Rin
        self.R2 = Rout
        self.height = height

    def depthLocal(self, r):
        r1 = r[:, 0::2]
        rn = norm(r1, axis=1)
        r1 = self.R1 - rn
        r2 = self.R2 - rn
        ins = np.array((r2 > 0) & (r1 < 0), dtype=int)
        return [r2, r1, ins]

    def cross(self, r, k):
        n = r.shape[0]
        [tin, tout] = tr.crossShellCyl(self.R1, self.R2, r, k)
        # apply other boundaries
        ty = tr.crossLayer(self.height, r, k, 1)
        for j in range(2):
            a = np.maximum(tin[:, j], ty[0])
            b = np.minimum(tout[:, j], ty[1])
            tin[:, j] = a[:, ]
            tout[:, j] = b[:, ]
        cm = []
        for i in range(n):
            tmp = [[tin[i, 0], tout[i, 0]], [tin[i, 1], tout[i, 1]]]
            cm.append(np.array(tmp))
        return cm

    def rayPaths(self, r, ki, kf):
        [ti1, ti2] = tr.crossShellCyl(self.R1, self.R2, r, ki)
        [tf1, tf2] = tr.crossShellCyl(self.R1, self.R2, r, kf)
        # apply other boundaries
        tiy = tr.crossLayer(self.height, r, ki, 1)
        for j in range(2):
            a = np.maximum(ti1[:, j], tiy[0])
            b = np.minimum(ti2[:, j], tiy[1])
            ti1[:, j] = a[:, ]
            ti2[:, j] = b[:, ]
        tfy = tr.crossLayer(self.height, r, kf, 1)
        for j in range(2):
            a = np.maximum(tf1[:, j], tfy[0])
            b = np.minimum(tf2[:, j], tfy[1])
            tf1[:, j] = a[:, ]
            tf2[:, j] = b[:, ]
        path1 = np.minimum(ti2, 0) - ti1
        path2 = tf2 - np.maximum(tf1, 0)
        mi = np.array(path1 > 0, dtype=int)
        mf = np.array(path2 > 0, dtype=int)
        res1 = np.sum(path1*mi, axis=1)*norm(ki)
        res2 = np.sum(path2*mf, axis=1)*norm(kf)
        return [res1, res2]

    def plotContours2(self, plt, proj, color, linestyle):
        theta = np.arange(0, 2*np.pi, 0.01)
        if ((proj == 0) or (proj == 2)):
            wx = self.R1
            hx = self.height*0.5
            x = [-wx, -wx, wx, wx, -wx]
            y = [-hx, hx, hx, -hx, -hx]
            plt.plot(x, y, color=color, linestyle=linestyle)
            wx = self.R2
            x = [-wx, -wx, wx, wx, -wx]
            plt.plot(x, y, color=color, linestyle=linestyle)
        elif (proj == 1):
            x1 = self.R1*np.cos(theta)
            y1 = self.R1*np.sin(theta)
            plt.plot(x1, y1, color=color, linestyle=linestyle)
            x2 = self.R2*np.cos(theta)
            y2 = self.R2*np.sin(theta)
            plt.plot(x2, y2, color=color, linestyle=linestyle)
            
    def plotContours(self, ax, proj, color, linestyle):
        theta = np.arange(-np.pi/2, np.pi/2, 0.01)
        gray = (0.2, 0.2, 0.2, 0.15)
        white = (1., 1., 1., 1.)
        if ((proj == 0) or (proj == 2)):
            w1 = self.R1
            w2 = self.R2
            h = self.height*0.5           
            # fill gray outer
            x = [-w2, -w2, w2, w2]
            y = [-h, h, h, -h]
            ax.fill_between([-w2, w2], [-h, -h], [h, h], facecolor=gray) 
            # fill white inner
            x = [-w1, -w1, w1, w1]
            y = [-h, h, h, -h]
            ax.fill_between([-w1, w1], [-h, -h], [h, h], facecolor=white)
            # contour outer rectangle
            x = [-w2, -w2, w2, w2, -w2]
            y = [-h, h, h, -h, -h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='a')
            # left line
            x = [-w1, -w1]
            y = [-h, h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='b')
            # right line
            x = [w1, w1]
            y = [-h, h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='c')
        elif (proj == 1):
            x = self.R2*np.sin(theta)
            # fill outer circle gray
            y2 = self.R2*np.cos(theta)
            y1 = -y2
            ax.fill_between(x, y1, y2, facecolor=gray)
            ax.plot(x, y1, color=color, linestyle=linestyle, label='a')
            ax.plot(x, y2, color=color, linestyle=linestyle, label='b')
            # fill inner circle white
            x = self.R1*np.sin(theta)
            y2 = self.R1*np.cos(theta)
            y1 = -y2
            ax.fill_between(x, y1, y2, facecolor=white)
            ax.plot(x, y1, color=color, linestyle=linestyle, label='c')
            ax.plot(x, y2, color=color, linestyle=linestyle, label='d')
