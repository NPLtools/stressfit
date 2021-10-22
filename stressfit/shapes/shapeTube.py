# -*- coding: utf-8 -*-
"""
Shape definition for a hollow cylinder
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapeTube(ShapeAbstract):
    """Define hollow cylinder with axis || y.
    
    The inner cyllindric hole can be off center, its position is 
    defined by ``ctr``. Define the reference surface by the value 
    of ``sref`` (0/1 for the inner/outer surface). The depth is
    calculated with respect to the reference cylindric surface 
    (regardless of the position along the axis). 

    Parameters
    ----------
    Rin : float
        Inner radius [mm].
    Rout : float
        Outer radius [mm].
    height : float
        Height in [mm].
    ctr : array_like(2) 
        Position (x,z) of the inner cyllinder centre.
    sref : int
        Index of the reference surface (0/1 for the inner/outer surface).

    """ 
    
    shape_type = 'cylinder_hollow'
    _srefs = ['inner', 'outer']
    def __init__(self, Rin=4.0, Rout=8.0, height=30.0, ctr=[0,0], sref=1):
        """Define hollow cylinder with axis || y."""
        super().__init__()
        self.R1 = min(Rin,Rout)
        self.R2 = max(Rin,Rout)
        self.R = [self.R1,self.R2] # for quick reference
        self.height = height
        i2 = [1, 0] # i2[self.sref] = index of the non-reference surface
        self.sref = sref
        self.uref = i2[self.sref]
        self._set_ctr(ctr)

    def _set_sref(self, sref):
        i2 = [1, 0] # i2[self.sref] = index of the non-reference surface
        old_sref = self.sref
        old_uref = self.uref
        if isinstance(sref, int):
            iref = sref
        else:
            iref = min(1,ShapeTube._srefs.index(sref))
        self.sref = max(0,min(1,iref))    
        self.uref = i2[self.sref]
        if old_sref != self.sref:
            sg = [-1, 1] # move relative to the reference circle
            ctr = sg[old_sref]*self.ctr[old_uref,:]
            self._set_ctr(ctr)

    def _set_ctr(self, ctr):
        """Set centre of the non-reference surface.
        
        The reference surface is considered centered in the local frame.
        """
        self.ctr = np.zeros((2,3))
        sg = [-1, 1] # move relative to the reference circle
        i2 = [1, 0] # i2[self.sref] = index of the non-reference surface
        j = i2[self.sref]
        self.ctr[j,:] = sg[self.sref]*np.array([ctr[0], 0, ctr[1]], 
                                               dtype=float) 

    def update(self, **kwargs):
        """Update parameters."""        
        if 'thickness' in kwargs:
            self.thickness = kwargs['thickness']
        if 'Rin' in kwargs:
            Rin = kwargs['Rin']
            Rout = self.R2
            self.R1 = min(Rin,Rout)
            self.R2 = max(Rin,Rout)
            self.R = [self.R1,self.R2] 
        if 'Rout' in kwargs:
            Rin = self.R1
            Rout = kwargs['Rout']
            self.R1 = min(Rin,Rout)
            self.R2 = max(Rin,Rout)
            self.R = [self.R1,self.R2] 
        if 'height' in kwargs:
            self.height = kwargs['height']
        if 'ctr' in kwargs:
            self._set_ctr(kwargs['ctr'])
        if 'sref' in kwargs:
            self._set_sref(kwargs['sref'])


    def depthLocal(self, r):
        """Calculate depths under the both cylindric surfaces.
        
        Reference surface is the one defined by the ``sref`` parameter
        (sref=0 denotes the inner surface). Calculates also a flag 
        marking points which are inside the object.
        
        **NOTE** Depth under the inner surface is measured as the radial 
        distance from the inner hole center, minus the inner radius. Hence
        the point is inside material if both depths are positive. 

        Parameters
        ----------
        r : array(:,3)
            Positions in local coordinates.

        Returns
        -------
        [depth1, depth2, inside]
            - depth1: depth under the reference surface
            - depth2: depth under the other surface
            - inside: (0|1) for outside|inside position
        """
        i2 = [1, 0] # index of the non-reference surface
        #r2 = r[:, 0::2] # get xz coordinates
        r20 = r[:, 0::2]-self.ctr[0,0::2]
        r21 = r[:, 0::2]-self.ctr[1,0::2]
        #sg = [-1, 1]
        # depths under the inner and outer surfaces
        depth = [norm(r20, axis=1)-self.R[0], self.R[1]-norm(r21, axis=1)]
        """
        for i in range(2):
            d = norm(r2-self.ctr[i,0::2], axis=1)
            
            depth[i] = sg[i]*(self.R[i] - d) 
        """
        # 'inside' flag
        qry = (abs(r[:,1])<0.5*self.height) & (depth[0] > 0) & (depth[1] > 0)
        ins = np.array(qry, dtype=int)
        return [depth[self.sref], depth[i2[self.sref]], ins]

    def cross(self, r, k):
        n = r.shape[0]
        #[tin, tout] = tr.crossShellCyl(self.R1, self.R2, r, k)
        [tin, tout] = tr.crossHollowCyl(self.R, self.ctr, r, k)
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
       # [ti1, ti2] = tr.crossShellCyl(self.R1, self.R2, r, ki)
       #[tf1, tf2] = tr.crossShellCyl(self.R1, self.R2, r, kf)
        [ti1, ti2] = tr.crossHollowCyl(self.R, self.ctr, r, ki)
        [tf1, tf2] = tr.crossHollowCyl(self.R, self.ctr, r, kf)
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
        # projection plane: 0: (z, y); 1: (x, y); 2: (x, z)
        if ((proj == 0) or (proj == 2)):
            if proj==0:
                dx1 = self.ctr[0,2]
                dx2 = self.ctr[1,2]
            else:
                dx1 = self.ctr[0,0]
                dx2 = self.ctr[1,0]
            w1 = self.R1
            w2 = self.R2
            h = self.height*0.5
            # fill gray outer
            ax.fill_between([-w2+dx2, w2+dx2], [-h, -h], [h, h], facecolor=gray) 
            # fill white inner
            ax.fill_between([-w1+dx1, w1+dx1], [-h, -h], [h, h], facecolor=white)
            # contour outer rectangle
            x = np.add([-w2, -w2, w2, w2, -w2],dx2)
            y = [-h, h, h, -h, -h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='a')
            # left line
            x = np.add([-w1, -w1], dx1)
            y = [-h, h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='b')
            # right line
            x = np.add([w1, w1], dx1)
            y = [-h, h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='c')
        elif (proj == 1):
            for i in [1,0]:
                self._fill_circle(ax, theta, i, color, linestyle)
            """
            # fill outer circle gray
            x = self.R2*np.sin(theta)+self.ctr[1,0]
            y2 = self.R2*np.cos(theta)+self.ctr[1,2]
            y1 = -y2+2*self.ctr[1,2]
            ax.fill_between(x, y1, y2, facecolor=gray)
            ax.plot(x, y1, color=color, linestyle=linestyle, label='a')
            ax.plot(x, y2, color=color, linestyle=linestyle, label='b')
            # fill inner circle white
            x = self.R1*np.sin(theta)+self.ctr[0,0]
            y2 = self.R1*np.cos(theta)+self.ctr[0,2]
            y1 = -y2+2*self.ctr[0,2]
            ax.fill_between(x, y1, y2, facecolor=white)
            ax.plot(x, y1, color=color, linestyle=linestyle, label='c')
            ax.plot(x, y2, color=color, linestyle=linestyle, label='d')
            """

    def _fill_circle(self, ax, theta, icirc, color, linestyle):
        gray = (0.2, 0.2, 0.2, 0.15)
        white = (1., 1., 1., 1.)
        colors = [white, gray]
        labels = [['a', 'b'],['c', 'd']]
        x = self.R[icirc]*np.sin(theta)+self.ctr[icirc,0]
        y2 = self.R[icirc]*np.cos(theta)+self.ctr[icirc,2]
        y1 = -y2+2*self.ctr[icirc,2]
        ax.fill_between(x, y1, y2, facecolor=colors[icirc])
        ax.plot(x, y1, color=color, linestyle=linestyle, label=labels[icirc][0])
        ax.plot(x, y2, color=color, linestyle=linestyle, label=labels[icirc][1])          
                  
        
        