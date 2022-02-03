# -*- coding: utf-8 -*-
"""
Extension of shapeTube. Defines multiple coaxial holes in a cylinder. 

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapeTubes(ShapeAbstract):
    """Define cylinder with axis || y and multiple coaxial holes.
    
    The inner cyllindric holes can be off center, their position is 
    defined by ``ctr``. The depth is
    calculated as the position along provided scan direction. 

    Parameters
    ----------
    Rout : float
        Outer radius [mm].
    Rin : list of float
        Inner radii [mm].
    height : float
        Height in [mm].
    ctr : list of array_like
        Positions (x,z) of the hole centres. One row per a hole. 
    sdir : array_like
        Scan direction in local coordinates.
    sctr : array_like
        Scan origin in local coordinates. 

    """ 
    
    shape_type = 'cylinder_holes'
    _srefs = ['inner', 'outer']
    def __init__(self, Rout=8.0, Rin=[4.0], height=30.0, ctr=[[0,0]], 
                 sdir=[0,0,1], sctr=[0,0,0]):
        """Define hollow cylinder with axis || y."""
        super().__init__()
        qry = np.all(np.array(Rin)<Rout)
        if not qry:
            raise Exception('Inner holes must be smaller than the outer radius.')
        self.R1 = Rin
        self.R2 = Rout
        self.R = np.array([self.R2]+self.R1) # for quick reference
        self.height = height
        self.set_scan(sdir, sctr)
        self.nh = len(Rin)
        qry = (self.nh==len(ctr))
        if not qry:
            raise Exception('Number of centers must equal the number of Rin values plus one.')
        self._set_ctr(ctr)
        
    def _set_ctr(self, ctr):
        """Set centre of the non-reference surface.
        
        The reference surface is considered centered in the local frame.
        """
        nc = len(ctr)
        self.ctr = np.zeros((nc+1,3))
        for i in range(1,nc+1):
            self.ctr[i,:] = np.array([ctr[i-1][0], 0, ctr[i-1][1]], dtype=float) 
            
    def update(self, **kwargs):
        """Update parameters."""        
        if 'Rin' in kwargs:
            self.R1 = kwargs['Rin']
            self.R = np.array([self.R2]+self.R1)
        if 'Rout' in kwargs:
            self.R2 = kwargs['Rout']
            self.R = np.array([self.R2]+self.R1) 
        if 'height' in kwargs:
            self.height = kwargs['height']
        if 'ctr' in kwargs:
            self._set_ctr(kwargs['ctr'])
        if 'sdir' in kwargs:
            self.set_scan(kwargs['sdir'], self._sctr)
        if 'sctr' in kwargs:
            self.set_scan(self._sdir, kwargs['sctr'])

    def depthLocal(self, r):
        """Calculate depths under the surfaces.
        
        Calculate depth as (i) projection on the scan vector and (ii) as
        depth under the surfaces of the sample and holes. The depth is 
        always positive towards the sample material, i.e. inside
        the sample and outside a hole.
        Calculates also a flag marking points which are inside the material
        (i.e. outside of al holes and inside the sample).
    
        Parameters
        ----------
        r : array(:,3)
            Positions in local coordinates.

        Returns
        -------
        [depth1, depth2, inside]
            - depth1: depth as projection on the scan vector
            - depth2: depths under the other surfaces
            - inside: (0|1) for outside|inside position
        """
        # projection of r on the scan vector
        #d = (r - self._sctr).dot(self._sdir)
        d = r.dot(self._sdir)
        nc = self.ctr.shape[0] 
        sg = np.ones(nc)
        sg[0] = -1
        depth = nc*[0]
        for i in range(nc):
            r2 = r[:,0::2] - self.ctr[i,0::2]
            depth[i] = sg[i]*(norm(r2, axis=1)-self.R[i])
        # 'inside' flag
        qry1 = np.all(np.array(depth) > 0, axis=0)
        
        qry = (abs(r[:,1])<0.5*self.height) & qry1
        ins = np.array(qry, dtype=int)
        return [d, depth, ins]

    def cross(self, r, k):
        n = r.shape[0]
        #[tin, tout] = tr.crossShellCyl(self.R1, self.R2, r, k)
        [tin, tout, ins] = tr.crossTubes(self.R, self.ctr, r, k)
        # apply other boundaries
        ty = tr.crossLayer(self.height, r, k, 1)
        for j in range(self.nh+1):
            a = np.maximum(tin[:, j], ty[0])
            b = np.minimum(tout[:, j], ty[1])
            tin[:, j] = a[:, ]
            tout[:, j] = b[:, ]
        cm = []
        for i in range(n):
            tmp = []
            for j in range(self.nh+1):
                tmp.append([tin[i, j], tout[i, j]])
            cm.append(np.array(tmp))
        return cm

    def rayPaths(self, r, ki, kf):
       # [ti1, ti2] = tr.crossShellCyl(self.R1, self.R2, r, ki)
       #[tf1, tf2] = tr.crossShellCyl(self.R1, self.R2, r, kf)
        [ti1, ti2, ins] = tr.crossTubes(self.R, self.ctr, r, ki)
        [tf1, tf2, ins] = tr.crossTubes(self.R, self.ctr, r, kf)
        # apply other boundaries
        tiy = tr.crossLayer(self.height, r, ki, 1)
        for j in range(self.nh+1):
            a = np.maximum(ti1[:, j], tiy[0])
            b = np.minimum(ti2[:, j], tiy[1])
            ti1[:, j] = a[:, ]
            ti2[:, j] = b[:, ]
        tfy = tr.crossLayer(self.height, r, kf, 1)
        for j in range(self.nh+1):
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


    def plotContours(self, ax, proj, color, linestyle):
        gray = (0.2, 0.2, 0.2, 0.15)
        white = (1., 1., 1., 1.)
        # projection plane: 0: (z, y); 1: (x, y); 2: (x, z)
        if ((proj == 0) or (proj == 2)):
            if proj==0:
                dx = self.ctr[:,2]
            else:
                dx = self.ctr[:,0]
            h = self.height*0.5
            w = self.R
            for i in range(len(self.R)):
                if i==0:
                    c = gray
                else:
                    c = white
                label = str(i)
                # fill 
                ax.fill_between([-w[i]+dx[i], w[i]+dx[i]], [-h, -h], [h, h], facecolor=c)
                # contour
                # left line
                x = np.add([-w[i], -w[i]], dx[i])
                y = [-h, h]
                ax.plot(x, y, color=color, linestyle=linestyle, label=label+'b')
                # right line
                x = np.add([w[i], w[i]], dx[i])
                ax.plot(x, y, color=color, linestyle=linestyle, label=label+'c')
            # outer contour
            x = np.add([-w[0], -w[0], w[0], w[0], -w[0]], dx[0])
            y = [-h, h, h, -h, -h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='a')
        elif (proj == 1):
            theta = np.arange(-np.pi/2, np.pi/2, 0.01)
            for i in range(len(self.R)):
                if i==0:
                    c = gray
                else:
                    c = white
                label = str(i)
                x = self.R[i]*np.sin(theta)+self.ctr[i,0]
                y2 = self.R[i]*np.cos(theta)+self.ctr[i,2]
                y1 = -y2+2*self.ctr[i,2]
                ax.fill_between(x, y1, y2, facecolor=c)
                ax.plot(x, y1, color=color, linestyle=linestyle, label='b'+label)
                ax.plot(x, y2, color=color, linestyle=linestyle, label='c'+label)

        