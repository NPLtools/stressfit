# -*- coding: utf-8 -*-
"""
Shape definition for a curved plate
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapePlateCurved(ShapeAbstract):
    """Define a bent plate.

    Plate with two independent, parabolic surfaces.
    Surface equation is z = 0.5*(rho[0]*x**2 + rho[1]*y**2) +- 0.5*thickness

    Parameters
    ----------
    thickness: float
        Plate z-thickness at (x,y) = 0
    length: float
        Plate x-length [mm]
    height: float
        Plate y-height [mm]
    rho1: array_like(2)
        Curvatures of the front surface, [hor, ver]
    rho2: array_like(2)
        Curvatures of the rear surface, [hor, ver]
    """
    shape_type = 'plate_curved'

    def __init__(self, thickness=5.0, length=50.0, height=15.0, 
                 rho1=[0.02, 0.0], rho2=[0.02, 0.0]):
        super().__init__()
        self.t = thickness
        self.rho1 = np.array(rho1)
        self.rho2 = np.array(rho2)
        self.size = np.array([length, height])

# overriden abstract methods:

    def update(self, **kwargs):
        """Update parameters."""
        if 'thickness' in kwargs:
            self.thickness = kwargs['thickness']
        if 'length' in kwargs:
            self.size = np.array([kwargs['length'],self.height])
        if 'height' in kwargs:
            self.size = np.array([self.length, kwargs['height']])
        if 'rho1' in kwargs:
            self.rho1 = np.array(kwargs['rho1'])
        if 'rho2' in kwargs:
            self.rho2 = np.array(kwargs['rho2'])
            
    def depthLocal(self, r):
        s1 = self.getSurface(r, 1)
        s2 = self.getSurface(r, -1)
        dep1 = s1 - r[:, 2]
        dep2 = s2 - r[:, 2]
        # 2nd order correction
        # NOTE: slows down execution by factor ~ 3
        N1 = self.getNormal(r, 1)
        N2 = self.getNormal(r, -1)
        dep1 *= N1[2, ]
        dep2 *= N2[2, ]
        lim = 0.5*self.size
        bounds1 = np.array(abs(r[:, 0]) < lim[0], dtype=int)
        bounds2 = np.array(abs(r[:, 1]) < lim[1], dtype=int)
        ins = np.array((dep1 > 0.) & (dep2 < 0.), dtype=int)
        return [dep1, dep2, ins*bounds1*bounds2]

    def cross(self, r, k):
        # get list of cross times [tin, tout]
        nr = r.shape[0]
        uno = np.ones(nr)
        li = tr.quad2D(1, -0.5*self.t, self.rho1, self.size, r, k)
        li2 = tr.quad2D(-1, 0.5*self.t, self.rho2, self.size, r, k)
        li.extend(li2)

        # collect times and signs into arrays at, asg
        i = 0
        t = []
        sg = []
        for X in li:
            if (np.mod(i, 2) == 0):
                t.append(li[i])
            else:
                sg.append(li[i])
            i += 1
        at = np.array(t).T
        asg = np.array(sg).T

        uno = np.ones(asg.shape)
        # get sorted entry times
        inf = 1.e12  # mask out other points by a large time value
        mskin = np.array((uno == asg), dtype=int)
        atin = mskin*at + (1-mskin)*inf
        atin.sort(axis=1)
        # get sorted exit times
        mskout = np.array((-uno == asg), dtype=int)
        atout = mskout*at + (1-mskout)*inf
        atout.sort(axis=1)

        # condition for a valid cross-point:
        # t_in < t_out, t < inf
        mask = np.array(
            (atout > atin) & (abs(atin) < inf) & (abs(atout) < inf), dtype=int)
        # valid entry and exit times, sorted
        tin = mask*atin - (1-mask)*inf
        tout = mask*atout + (1-mask)*inf

        # apply other boundaries
        tx = tr.crossLayer(self.size[0], r, k, 0)
        ty = tr.crossLayer(self.size[1], r, k, 1)
        tmi = np.minimum(tx[1], ty[1])
        tma = np.maximum(tx[0], ty[0])
        for j in range(tin.shape[1]):
            tin[:, j] = np.maximum(tin[:, j], tma.T)
        for j in range(tout.shape[1]):
            tout[:, j] = np.minimum(tout[:, j], tmi.T)

        # get list of arrays with cross-section times
        # one array per position
        # each array contains entry and exit times on rows
        nt = tin.shape
        cm = []
        for j in range(nt[0]):
            tmp = []
            for m in range(nt[1]):
                if (mask[j, m] and (tin[j, m] < tout[j, m])):
                    tmp.append([tin[j, m], tout[j, m]])
            cm.append(np.array(tmp))
        return cm

    def rayPaths(self, r, ki, kf):
        n = r.shape[0]
        ski = np.sign(ki[2])
        skf = np.sign(kf[2])
        di = (-ski*0.5*self.t - r[:, 2]).reshape((n, 1))
        df = (skf*0.5*self.t - r[:, 2]).reshape((n, 1))
        ti = di/ki[2]
        tf = df/kf[2]
        # now we should have ti<0 and tf>0
        # TODO: exclude exceptions

        # r + t*k = position at flat boundaries, z = +- t/2
        # 1st order correction for bending,
        ri = r + ti*ki
        rf = r + tf*kf
        # distance to the surface, corrected
        di = self.getSurface(ri, -ski) - r[:, 2]
        df = self.getSurface(rf, skf) - r[:, 2]

        # time to surface
        ti = di/ki[2]
        tf = df/kf[2]

        # clip on size
        ski = np.sign(ki[0:2])
        skf = np.sign(kf[0:2])
        for i in range(2):
            dx = -ski[i]*0.5*self.size[i] - r[:, i]
            if ski[i] != 0:
                ti = np.maximum(dx/ki[i], ti)
            dx = skf[i]*0.5*self.size[i] - r[:, i]
            if skf[i] != 0:
                tf = np.minimum(dx/kf[i], tf)
        # paths to surface
        sg1 = np.array(ti < 0, dtype=int)
        sg2 = np.array(tf > 0, dtype=int)
        path1 = -sg1*ti*norm(ki)
        path2 = sg2*tf*norm(kf)
        return [path1, path2]

# class specific methods:

    def getSurface(self, r, isurf):
        x2 = r[:, 0]**2
        y2 = r[:, 1]**2
        s = 0.5*(1+np.sign(isurf))
        rh = s*self.rho2+(1-s)*self.rho1
        y = 0.5*(rh[0]*x2 + rh[1]*y2)
        y += np.sign(isurf)*0.5*self.t
        return y

    def getNormal(self, r, isurf):
        uno = np.ones(r.shape[0])
        # d = self.depthLocal(r)
        s = 0.5*(1+np.sign(isurf))
        rh = s*self.rho2 + (1-s)*self.rho1
        v = np.array([-rh[0]*r[:, 0], -rh[1]*r[:, 1], uno])
        N = v/norm(v, axis=0)
        return N

    def getContours(self,rang, n):
        y = np.zeros((n,2))
        x = np.linspace(rang[0], rang[1], num=n, endpoint=True)
        yy = np.zeros(n)
        r = np.array([x, yy]).reshape((n,2))
        y[:,0]=self.getSurface(r, -1)
        y[:,1]=self.getSurface(r, 1)
        return [x, y]    

    def plotContours2(self, plt, proj, color, linestyle):
        ax = plt.axes()
        xlim=ax.get_xlim()
        n = 200
        y = np.zeros((n,2))
        x = np.linspace(xlim[0], xlim[1], num=n, endpoint=True)
        yy = np.zeros(n)
        r = np.array([x, yy]).reshape((n,2))
        y[:,0] =  self.getSurface(r, -1)
        y[:,1] =  self.getSurface(r, 1)
        for i in range(2):
            plt.plot(x, y[:,i],color=color, linestyle=linestyle, marker='none')   

    def plotContours(self, ax, proj, color, linestyle):
        gray = (0.2, 0.2, 0.2, 0.15)
        n = 100
        # z,y
        if (proj == 0):
            y = np.linspace(-0.5*self.size[1], 0.5*self.size[1], num=n, endpoint=True)
            xx = np.zeros(n)
            r = np.array([xx, y]).T
            z = np.zeros((n,2))
            z[:,0] =  self.getSurface(r, -1)
            z[:,1] =  self.getSurface(r, 1)
            # fill gray
            ax.fill_betweenx(y, z[:,0], x2=z[:,1], facecolor=gray) 
            ax.plot(z[:,0], y,
                    color=color, linestyle=linestyle, marker=None, label='c')
            ax.plot(z[:,1], y,
                    color=color, linestyle=linestyle, marker=None, label='d')
            ax.plot([z[0,0],z[0,1]], [-0.5*self.size[1],-0.5*self.size[1]],
                    color=color, linestyle=linestyle, marker=None, label='a')
            ax.plot([z[-1,0],z[-1,1]], [0.5*self.size[1],0.5*self.size[1]],
                    color=color, linestyle=linestyle, marker=None, label='b')       
        # x,z
        elif (proj == 1):
            x = np.linspace(-0.5*self.size[0], 0.5*self.size[0], num=n, endpoint=True)
            yy = np.zeros(n)
            r = np.array([x, yy]).T
            z = np.zeros((n,2))
            z[:,0] =  self.getSurface(r, -1)
            z[:,1] =  self.getSurface(r, 1)
            # fill gray
            ax.fill_between(x, z[:,0], y2=z[:,1], facecolor=gray) 
            ax.plot(x, z[:,0],
                    color=color, linestyle=linestyle, marker=None, label='c')
            ax.plot(x, z[:,1],
                    color=color, linestyle=linestyle, marker=None, label='d')
            ax.plot([-0.5*self.size[0], -0.5*self.size[0]], [z[0,0],z[0,1]], 
                    color=color, linestyle=linestyle, marker=None, label='a')
            ax.plot([0.5*self.size[0], 0.5*self.size[0]], [z[-1,0],z[-1,1]],
                    color=color, linestyle=linestyle, marker=None, label='b')

        # x,y
        elif (proj == 2):
            x1 = np.array([0.5*self.size[0], 0.5*self.size[0]])
            x = np.array([-0.5*self.size[0], 0.5*self.size[0]])
            y1 = np.array([-0.5*self.size[1], -0.5*self.size[1]])
            y2 = np.array([0.5*self.size[1], 0.5*self.size[1]])
            # fill gray
            ax.fill_between(x, y1, y2=y2, facecolor=gray) 
            ax.plot(x, y1, color=color, linestyle=linestyle, marker=None, label='c')
            ax.plot(x, y2, color=color, linestyle=linestyle, marker=None, label='d')
            ax.plot(-x1, y1, color=color, linestyle=linestyle, marker=None, label='a')
            ax.plot(x1, y2, color=color, linestyle=linestyle, marker=None, label='b')
