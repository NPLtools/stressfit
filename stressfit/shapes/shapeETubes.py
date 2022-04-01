# -*- coding: utf-8 -*-
"""
Extension of shapeTube. Defines multiple coaxial holes in a cylinder, with 
elliptic basis. 

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from numpy.linalg import norm
from matplotlib.patches import Ellipse
from . import tracing as tr
from .shapeAbstract import ShapeAbstract


class ShapeETubes(ShapeAbstract):
    """Define cylinder with axis || y and multiple coaxial elliptic holes.
    
    The inner holes can be at any position, with elliptic basis of 
    any orientation. 
    The depth is calculated as the position projected on the scan direction. 
    The configuration of holes is defined as a list of dict parameters.
    
    The holes should not overlap with each other or with the sample surface.
    No verification of this constraint is made!

    Parameters
    ----------
    a : float
        First semi-axis of the sample basis, parallel to x if angle=0. 
    b : float
        Second semi-axis of the sample basis, parallel to z if angle=0.  
    angle : float
        Angle of the semiaxis a with respect to x [deg]
    height : float
        Height in [mm].
    holes : list of dict
        Each item defines a hole. 
        The parameters are
            x, z : coordinates of the centre
            a, b : ellipse semi-axes
            angle : angle of the semiaxis a with respect to x [deg].
    sdir : array_like
        Scan direction in local coordinates.
    sctr : array_like
        Scan origin in local coordinates. 
    """ 
    
    shape_type = 'cylinder_eholes'
    def __init__(self, a=8.0, b=8.0, angle=0.0, height=30.0, holes=[], 
                 sdir=[0,0,1], sctr=[0,0,0]):
        """Define hollow cylinder with axis || y."""
        super().__init__()       
        self.a = a
        self.b = b
        self.angle = angle*np.pi/180    
        self.height = height
        self.holes = holes
        self._set_holes(holes)
        self.set_scan(sdir, sctr) 
        
    def get_param(self):
        out = {}
        out['a'] = self.a
        out['b'] = self.b
        out['height'] = self.height
        out['angle'] = self.angle*180/np.pi
        out['holes'] = self.holes
        return out
    
    def _set_holes(self, holes):
        """Set parameters of holes.
        
        Parameters
        ----------
        holes : list of dict
            Each item defines a hole. 
            The parameters are:
                - x, z : coordinates of the centre 
                - a, b : ellipse semi-axes
                - angle : angle of the semi-axis a with respect to x [deg].
        """
        if self.holes != holes:
            self.holes = holes.copy()
        self.nh = len(self.holes)
        self.ctr = np.zeros((self.nh+1,3))
        self.gm = (self.nh+1)*[0]
        self.gm[0] = tr.get_metric_ell2D(self.a,self.b,self.angle) 
        for i in range(self.nh):
            h = holes[i]
            self.ctr[i+1,:] = np.array([h['x'], 0., h['z']], dtype=float)
            self.gm[i+1] = tr.get_metric_ell2D(h['a'],h['b'],h['angle']*np.pi/180) 
                    
    def update(self, **kwargs):
        """Update parameters."""    
        _calc = False
        if 'a' in kwargs:
            self.a = kwargs['a']
            _calc = True
        if 'b' in kwargs:
            self.b = kwargs['b']
            _calc = True
        if 'angle' in kwargs:
            self.angle = kwargs['angle']*np.pi/180
            _calc = True
        if 'height' in kwargs:
            self.height = kwargs['height']
        if 'sdir' in kwargs:
            self.set_scan(kwargs['sdir'], self._sctr)
        if 'sctr' in kwargs:
            self.set_scan(self._sdir, kwargs['sctr'])
        if 'holes' in kwargs:
            self._set_holes(kwargs['holes'])
        elif _calc:
            self._set_holes(self.holes)                 

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
            rsq = np.sum(r2*(r2.dot(self.gm[i])), axis=1)
            rnorm = norm(r2, axis=1)
            depth[i] = sg[i]*(np.sqrt(rsq)-1.0)*rnorm
        # 'inside' flag
        qry1 = np.all(np.array(depth) > 0, axis=0)
        qry = (abs(r[:,1])<0.5*self.height) & qry1
        ins = np.array(qry, dtype=int)
        return [d, depth, ins]

    def cross(self, r, k):
        n = r.shape[0]
        [tin, tout, ins] = tr.crossETubes(self.gm, self.ctr, r, k)
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
        [ti1, ti2, ins] = tr.crossETubes(self.gm, self.ctr, r, ki)
        [tf1, tf2, ins] = tr.crossETubes(self.gm, self.ctr, r, kf)
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
        _deg = np.pi/180
        
        def par(i):
            if i==0:
                c = gray
                a = self.a
                b = self.b
                angle = self.angle
            else:
                c = white
                hole =  self.holes[i-1]
                a = hole['a']
                b = hole['b']
                angle = hole['angle']*_deg
            return [c, a, b, angle]
        
        # projection plane: 0: (z, y); 1: (x, y); 2: (x, z)
        if ((proj == 0) or (proj == 2)):
            if proj==0:
                dx = self.ctr[:,2]
                w0 = self.b
            else:
                dx = self.ctr[:,0]
                w0 = self.a
            h = self.height*0.5
            for i in range(len(self.gm)):
                [c, a, b, angle] = par(i)
                label = str(i)
                [wx, wy] = tr.ell_limits(a, b, angle)
                if proj==0:
                    w = wy
                else:
                    w = wx
                # fill 
                ax.fill_between([-w+dx[i], w+dx[i]], [-h, -h], [h, h], facecolor=c)
                # contour
                # left line
                x = np.add([-w, -w], dx[i])
                y = [-h, h]
                ax.plot(x, y, color=color, linestyle=linestyle, label=label+'b')
                # right line
                x = np.add([w, w], dx[i])
                ax.plot(x, y, color=color, linestyle=linestyle, label=label+'c')
            # outer contour
            x = np.add([-w0, -w0, w0, w0, -w0], dx[0])
            y = [-h, h, h, -h, -h]
            ax.plot(x, y, color=color, linestyle=linestyle, label='a')
        elif (proj == 1):
            for i in range(len(self.gm)):
                [c, a, b, angle] = par(i)
                label = str(i)
                ctr = self.ctr[i,0::2]
                # Ellipse does not work, it hides scatter plot ...
                #ell = Ellipse(xy=ctr, width=2*a, height=2*b, angle=-angle/_deg, 
                #              edgecolor=color, linestyle=linestyle, facecolor=c)
                #ax.add_patch(ell)
                _plot_ellipse(ax, ctr=ctr, a=a, b=b, angle=angle, 
                              edgecolor=color, facecolor=c, 
                              linestyle=linestyle, label=str(i))
    
def _plot_ellipse(ax, ctr=[0,0], a=1, b=1, angle=0, edgecolor='b', 
                  facecolor='white', linestyle='-', label=''):
    thx = np.arctan2(b*np.sin(angle),a*np.cos(angle))
    theta = np.arange(thx, thx+np.pi, 0.01)
    [x,y1] = tr.ell(theta, a, b, angle)
    y2 = -y1[::-1]
    ax.fill_between(x+ctr[0], y1+ctr[1], y2+ctr[1], facecolor=facecolor)
    ax.plot(x+ctr[0], y1+ctr[1], color=edgecolor, linestyle=linestyle, 
            label='b'+label)
    ax.plot(x+ctr[0], y2+ctr[1], color=edgecolor, linestyle=linestyle, 
            label='c'+label)
        