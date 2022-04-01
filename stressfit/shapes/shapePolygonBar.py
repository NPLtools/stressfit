# -*- coding: utf-8 -*-
"""
Shape definition for a column with polygonal basis.

@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from numpy.linalg import norm
from . import tracing as tr
from .shapeAbstract import ShapeAbstract

_eps = 1.0e-12


class ShapePolygonBar(ShapeAbstract):
    """Define a bar with polygonal basis (axis along y).
    
    Define the basis as a regular polygon by giving the number of edges. 
    At zero position, the centre of the polygon is at the origin of local 
    coordinates.
    
    Use the optional `edges` parameter to define any (irregular) polygon 
    instead.
    
    Parameters
    ----------
    num: int
        Plate thickness in [mm].
    side: float
        Length of one side [mm].
    height : float
        Height along y [mm].
    angle: float
        Rotation of the basis [deg]. Zero rotation sets the one side
        parallel to x-axis.
    edges : list
        Optional, allows to set x,z coordinates of the edges individually.
        If defined, it overrides num, side and angle parameters. 
    sdir : array_like
        Scan direction in local coordinates.
    sctr : array_like
        Scan origin in local coordinates.         
        
    """
    shape_type = 'polygon_bar'

    def __init__(self, num=6, side=10.0, height=50.0, angle=0.0, edges=None,
                 sdir=[0,0,1], sctr=[0,0,0]):
        super().__init__()
        self.num = int(num)
        self.side = side
        self.height = height
        self.angle = angle
        self.set_scan(sdir, sctr)
        self._define_edges(edges=edges)
    
        
    def _define_edges(self, edges=None):
        """Calculate edges coordinates from input parameters."""
        self._input_edges = edges
        if edges is not None:
            self.num = len(edges)
            self.angle = 0.0
            self.edges = []
            for i in range(self.num):
                self.edges.append([edges[i][0], 0, edges[i][1]])
            self._def_by_edges = True
        else:
            deg = np.pi/180
            a = 360/self.num
            R = 0.5*self.side/np.sin(0.5*a*deg)
            self.edges = []
            for i in range(self.num):
                ang = (i+0.5)*a + self.angle
                x = R*np.sin(ang*deg)
                z = R*np.cos(ang*deg)
                self.edges.append([x, 0, z])
            self._def_by_edges = None
        # calculate sides, centres and normals and side axes of the walls 
        self.sides = []
        self.centres = []
        self.normals = []
        self.axes = []
        for i in range(len(self.edges)):
            # vector along side
            rs = np.subtract(self.edges[i],self.edges[i-1])
            slen = np.sqrt(rs[0]**2+rs[2]**2)
            rs = rs/slen
            ry = np.array([0.,1.,0.])
            # unit vector normal to the side/wall (towards inside)
            rn = np.cross(ry,rs)
            wnorm = rn/np.sqrt(rn[0]**2+rn[2]**2)
            # wall centre
            wctr =  0.5*np.add(self.edges[i],self.edges[i-1])
            # side length
            self.sides.append(slen)
            self.centres.append(wctr)
            self.normals.append(wnorm)
            self.axes.append(rs)
    
    def get_param(self):
        out = {'height': self.height}
        if self._def_by_edges:
            out['edges'] = self._input_edges
        else:
            out['num'] = int(self.num)
            out['side'] = self.side
            out['angle'] = self.angle
        return out

    def update(self, **kwargs):
        """Update parameters."""
        if 'num' in kwargs:
            self.num = int(kwargs['num'])
        if 'side' in kwargs:
            self.side = kwargs['side']
        if 'height' in kwargs:
            self.height = kwargs['height']
        if 'angle' in kwargs:
            self.angle = kwargs['angle']
        if 'sdir' in kwargs:
            self.set_scan(kwargs['sdir'], self._sctr)
        if 'sctr' in kwargs:
            self.set_scan(self._sdir, kwargs['sctr'])            
        if 'edges' in kwargs:
            self._define_edges(edges=kwargs['edges'])
        else:
            self._define_edges()
            
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
        d = r.dot(self._sdir)
        
        nw = self.num
        d2 = (nw+2)*[0]
        # normal distances from wall centres
        for i in range(nw):
            d2[i] = (r - self.centres[i]).dot(self.normals[i])
        # add distances from the top and bottom wall
        rh = np.array([0.0, 0.5*self.height, 0.0])
        kh = np.array([0.0, -1.0, 0.0])
        d2[nw] = (r - rh).dot(kh)
        d2[nw+1] = (r + rh).dot(-kh)
        # is inside if all distances are positive
        q = np.array(d2) > 0
        ins = np.array(np.all(q, axis=0), dtype=int)
        return [d, d2, ins]

    def cross(self, r, k):
        t = tr.crossPolygon(r, k, self.height, self.sides, self.centres, 
                         self.normals, self.axes)
        tmp = np.array(t).T
        return list(tmp)

    def rayPaths(self, r, ki, kf):
        [tin, ti2] = tr.crossPolygon(r, ki, self.height, self.sides, self.centres, 
                         self.normals, self.axes)
        [tf1, tout] = tr.crossPolygon(r, kf, self.height, self.sides, self.centres, 
                         self.normals, self.axes)
        
        sg1 = np.array((tin < 0) & (tin > -tr._inf), dtype=int)
        sg2 = np.array((tout > 0) & (tout < tr._inf), dtype=int)
        path1 = -sg1*tin*norm(ki)
        path2 = sg2*tout*norm(kf)
        return [path1, path2]

    def plotContours(self, ax, proj, color, linestyle):
         # projection plane: 0: (z, y); 1: (x, z); 2: (x, y)
        _eps = 1e-6
        gray = (0.2, 0.2, 0.2, 0.15)
        if ((proj == 0) or (proj == 2)):
        # side view
            if proj==0:
                ix = 2
                iy = 0
                sg = 1
            else:
                ix = 0
                iy = 2
                sg = -1
            xl = []
            yl = []
            # collect ix coordinates for visible edges 
            for i in range(self.num):
                e = self.edges[i]
                xl.append(e[ix])
                yl.append(e[iy])
            xmi = min(xl)
            xma = max(xl)
            lines = []
            for i in range(len(xl)):
                if sg*yl[i]>-_eps or xl[i] == xmi or xl[i] == xma:
                    lines.append(xl[i])
            h = 0.5*self.height
            # fill gray the visible surface
            ax.fill_between([xmi, xma], [-h, -h], [h, h], facecolor=gray) 
            # draw vertical lines = edges
            for i in range(len(lines)):
                x = lines[i]
                ax.plot([x, x], [-h, h], color=color, linestyle=linestyle, 
                        marker=None, label='a'+str(i))
            # top/bottom contour
            ax.plot([xmi, xma], [-h, -h], color=color, linestyle=linestyle, 
                    marker=None, label='bottom')
            ax.plot([xmi, xma], [h, h], color=color, linestyle=linestyle, 
                    marker=None, label='top')
        else:
        # top view - draw the polygon
            x = []
            y = []
            for e in self.edges:
                x.append(e[0])
                y.append(e[2])
            # close polygon
            x.append(x[0]) 
            y.append(y[0]) 
            # fill gray
            ax.fill(x, y, facecolor=gray)
            # plot contour
            ax.plot(x, y, color=color, linestyle=linestyle, marker=None, 
                    label='a')



#%%


