# -*- coding: utf-8 -*-
"""
Handles scan geometry and MC integration methods.

Created on Sat Sep 16 15:36:02 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from .cuda import cp, asnumpy
from .primitives import Transform

# TODO introduce initial transform for scan origin - simplify positioning

class ScanGeometry():
    """Scan geometry class.
    
    Defines linear and circular scans with methods for coordinate conversions.
    
    Linear scan:
        Define an array of translation steps, unit step size and origin.
        Use constructor :meth:`ScanGeometry.linear`. Step size is in [mm].
    Circular scan:
        Define an array of angular steps, unit step size, rotation origin and
        rotation axis. Step size is in [deg]. 
        Use constructor :meth:`ScanGeometry.circular`.
    
    Note that any step sequence can be defined by using non-integer values in 
    the steps array.
    
    Example
    -------
    .. highlight:: python
    .. code-block:: python
    
        # Define a linear scan with three positions at x = [-0.9, 0.1, 1.7], 
        # y=z=0.
        scan = ScanGeometry.linear(steps=[-0.5, 0, 0.8], step_size=2, 
                                   origin=[0.1, 0, 0])

        # Define a scan by rotation of y-axis from -60 to 60 deg, rotation 
        # centre at [0,0,-10] and origin at [0,0,5]
        steps = np.linspace(-60, 60, num=13)
        sscan = ScanGeometry.circular(steps, origin=[0,0,5], axis=[0,1,0], 
                                      origin_rot=[0,0,-10])

    """
    
    LINEAR = 0
    CIRCULAR = 1
    
    def __init__(self, shape:int, steps, origin=[0,0,0], axis=[1,0,0],
                 step_size=1.0, origin_rot=[0,0,0]):
        if not shape in [0,1]:
            raise Exception('Unsupported scan shape: {}'.format(shape))
        self._shape = shape
        self._steps = cp.array(steps, dtype=float)
        self._origin = cp.array(origin, dtype=float)
        self._step_size = step_size
        a = cp.array(axis, dtype=float)
        self._axis = a/cp.linalg.norm(a) # normalize direction axis
        self._origin_rot = cp.array(origin_rot, dtype=float)

    @classmethod
    def linear(cls, steps, origin=[0,0,0], axis=[0,1,0], step_size=1.0):
        """Define a linear scan.
        
        Parameters
        ----------
        steps : array_like
            Step positions in [mm]
        origin : array_like
            Coordinates of the scan origin (step=0)
        axis : array_like
            Scan direction
        step_size : float
            Step size factor (multiplies the values of `steps`)
        
        """
        c = cls(ScanGeometry.LINEAR, steps, origin=origin, axis=axis, 
                step_size=step_size)
        return c

    @classmethod
    def circular(cls, steps, origin=[1,0,0], axis=[0,1,0], origin_rot=[0,0,0],
                 step_size=1.0):
        """Define a circular scan.
        
        Parameters
        ----------
        steps : array_like
            Step positions in [deg]
        origin : array_like
            Coordinates of the scan origin (step=0)
        axis : array_like
            Rotation axis direction
        origin_rot : array_like
            Origin of the rotation axis (a point on the axis)
        step_size : float
            Step size factor (multiplies the values of `steps`)
        
        """
        c = cls(ScanGeometry.CIRCULAR, steps, origin=origin, axis=axis, 
                step_size=step_size, origin_rot=origin_rot)
        return c

    @property
    def nsteps(self):
        """Number of scan steps."""
        return len(self._steps)

    @property
    def shape(self):
        """Scan shape."""
        return self._shape

    @property
    def steps(self):
        """Step array."""
        return self._steps

    @property
    def units(self):
        """Scan step units."""
        if self.shape == ScanGeometry.CIRCULAR:
            return 'deg'
        else:
            return 'mm'

    def rot_axes(self):
        """Return carthesian basis vectors for rotation scan.
        
        Returns
        -------
        x : array_like
            Vector connecting rotation axis and scan origin 
            (radial component).
        y : array_like
            Rotation axis
        z : array_like
            Normal to x and y (hoop component)
        """
        y = self._axis # axis is already normalized
        z1 = cp.cross(self._origin - self._origin_rot,y)
        x1 = cp.cross(y,z1)
        x = x1/cp.linalg.norm(x1)
        z = z1/cp.linalg.norm(z1)
        return x,y,z

    def lin_axes(self):
        """Return carthesian basis vectors for linear scan.
        
        Returns
        -------
        x : array_like
            Scan direction.
        y : array_like
            normal to x and z
        z : array_like
            Normal to scan direction and [0,1,0]. 
            If scan direction is near [0,1,0], then normal to [1,0,0].
        """
        x = self._axis 
        y = cp.array([0,1,0], dtype=float)
        z1 = cp.cross(x,y)
        z0 = cp.linalg.norm(z1)
        if abs(z0)<1e-3:
            y = cp.array([1,0,0], dtype=float)
            z1 = cp.cross(x,y)
            z0 = cp.linalg.norm(z1)
        z = z1/z0
        y = cp.cross(z,x)
        return x,y,z

    def axes(self):
        """Return carthesian basis for scan coordinates.
        
        Calls :meth:`lin_axes` and :meth:`rot_axes` for linear and rotation
        scans, respectively.
        """
        if self._shape == ScanGeometry.CIRCULAR:
            x,y,z = self.rot_axes()
        else:
            x,y,z = self.lin_axes()
        return x,y,z

    def _positions_rot(self):
        deg = cp.pi/180
        x,y,z = self.rot_axes()
        rmat = cp.array([x, y, z])       
        pts = []
        trans = []
        r0 = cp.zeros(3)
        for i in range(len(self._steps)):
            s = self._step_size*self._steps[i]*deg
            co = cp.cos(s)
            si = cp.sin(s)
            rm = cp.eye(3)
            rm[0,2] = -si
            rm[2,0] = si
            rm[0,0] = co
            rm[2,2] = co
            rmat1 = rmat.T.dot(rm)
            R = rmat1.dot(rmat)
            tr = Transform(rmat=R, orig=-self._origin,
                           rotctr=self._origin_rot, order=1)
           # o = rmat1.dot(r1) + self._origin_rot
            o = tr.r_to_loc(r0)
            pts.append(o)
            trans.append(tr)
        return {'r': cp.array(pts).T, 'tr':trans}

    def _positions_lin(self):
        r = self._step_size*cp.outer(self._steps,self._axis) + self._origin
        trans = []        
        for i in range(len(self._steps)):
            tr = Transform(orig=-r[i], order=1)
            trans.append(tr)
        return {'r': r.T, 'tr':trans}

    def positions(self):
        """Calculate transformations which generate scan positions.
        
        Returns
        -------
        dict
            r : array_like
                Translation vectors corresponding to scan positions,
                shape = (3,:).
            tr : list
                A list of Transform objects, which generate positions and 
                direction vectors at each scan step.         
        
        Usage
        -----
        .. highlight:: python
        .. code-block:: python
        
            P = scan.positions() # calculate scan positions 
            
            # get position of the i-th scan point
            r = P['r'][:,i]
            
            # alternatively, calculate it using the Transform object
            r0 = cp.zeros(3)
            r = P['tr'][i].r_to_loc(r0)
            
            # transform velocity vector v0 from scan origin to the i-th scan point
            v = P['tr'][i].v_to_loc(v0)
            
            # For linear scans, it is simply:
            v = v0
            r = P['tr'][i].iorig        

        """    
        res = None
        if self._shape == ScanGeometry.LINEAR:
            res = self._positions_lin()
        elif self._shape == ScanGeometry.CIRCULAR:
            res = self._positions_rot()
        return res


    def to_scan_coord_rot(self, r):
        """Convert r to rotation scan cordinates.
        
        Uses :meth:`rot_axes` to define reference carthesian basis.
        
        Returns r in cylindric coordinates: r (radial), theta (angle [deg])
        and h (axial)
        
        Returns
        -------
        dict
            r, theta, h, pos : 
                r converted to circular scan coordinates,
                pos = theta, position along scan.
        
        """
        if self._shape != ScanGeometry.CIRCULAR:
            return None
        x,y,z = self.rot_axes()        
        rc = cp.asarray(r) - self._origin_rot.reshape(3,1)
        rx = cp.dot(x,rc)
        ry = cp.dot(y,rc)
        rz = cp.dot(z,rc)
        r0 = cp.sqrt(rx**2+rz**2)
        a = cp.arctan2(rx, rz)
        o = self._origin - self._origin_rot        
        a0 = cp.arctan2(cp.dot(x,o), cp.dot(z,o))
        theta = a - a0
        # put theta on -pi,+pi interval
        q1 = cp.array(theta<-cp.pi, dtype=int)
        q2 = cp.array(theta>cp.pi, dtype=int)
        theta = theta + q1*2*cp.pi
        theta = theta - q2*2*cp.pi
        res = {}
        res['r'] = r0
        res['theta'] = theta*180/cp.pi
        res['h'] = ry
        res['pos'] = res['theta']
        return res

    def to_scan_coord_lin(self, r):
        """Convert r to linear scan cordinates.
        
        Uses :meth:`lin_axes` to define reference carthesian basis.
        
        Returns
        -------
        dict
            x, y, z, pos : 
                r converted to linear scan coordinates,
                pos = x, position along scan.
        
        """
        if self._shape != ScanGeometry.LINEAR:
            return None
        x,y,z = self.lin_axes()
        rx = cp.dot(x,r)
        ry = cp.dot(y,r)
        rz = cp.dot(z,r)
        x0 = cp.dot(self._origin,r)
        res = {}
        res['x'] = rx - x0
        res['y'] = ry
        res['z'] = rz
        res['pos'] = res['x']
        return res

    def to_scan_coord(self, r):
        """Convert r to scan cordinates.
        
        Calls :meth:`to_scan_coord_lin` and :meth:`to_scan_coord_rot` for 
        linear and rotation scans, respectively.
        
        Returns
        -------
        dict
            pos : array_like
                Position along scan.
        
        Other returned items depend on the scan shape.
        
        """
        if self._shape == ScanGeometry.CIRCULAR:
            return self.to_scan_coord_rot(r)
        else:
            return self.to_scan_coord_lin(r)


#%% Tests
# TODO extract tests to separate module
def test_scan():
    """Test circular and linear scan settings."""
    from matplotlib import pyplot as plt
    _eps = 1e-15
# test circular scan
    step_size = 1
    steps = cp.linspace(0,120,num=9)
    origin = [1,0,0]
    axis = [0,0,1]
    origin_rot = [-1,0,0]
    print('Testing circular scan ... '.format(),end='')
    scan = ScanGeometry.circular(steps,origin=origin, axis=axis, 
                                 origin_rot=origin_rot, step_size=step_size)
    pos = scan.positions()
    # test points
    r = asnumpy(pos['r'])
    # attached vectors along x and y
    x_ax = np.array([1,0,0])
    z_ax = np.array([0,1,0])
    a = {'x':[], 'z':[]}
    for i in range(scan.nsteps):
        x_axt = asnumpy(pos['tr'][i].v_to_loc(x_ax))
        z_axt = asnumpy(pos['tr'][i].v_to_loc(z_ax))
        vx = np.array([r[:,i], r[:,i]+x_axt]).T
        vz = np.array([r[:,i], r[:,i]+z_axt]).T
        a['x'].append(vx)
        a['z'].append(vz)
    # plot scan points with attached vectors
    xlim = [-3,3]
    ylim = [-3,3]
    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.plot(r[0], r[1], color='black', marker='x', linestyle='none', markersize=3)    
    ax.plot(origin[0], origin[1], color='red', marker='o', linestyle='none', 
            markersize=5)
    ax.plot(origin_rot[0], origin_rot[1], color='red', marker='^', 
            linestyle='none', markersize=5)
    for i in range(scan.nsteps):
        ax.plot(a['x'][i][0], a['x'][i][1], color='red', linestyle='-', 
                linewidth=0.5)
        ax.plot(a['z'][i][0], a['z'][i][1], color='blue', linestyle='-', 
                linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid()
    title = 'origin: '+3*'{:g} '+'axis: '+3*'{:g} '+'step: {:g}'
    plt.title(title.format(*origin, *axis, step_size), fontsize = 8)    
    ax.invert_xaxis()
    plt.show()     
    # assert value
    r0 = np.array([0,0,0])
    r_test = pos['tr'][6].r_to_loc(r0)
    rdiff = asnumpy(r_test) - np.array([-1,2,0])
    qry = np.sum(rdiff**2)
    assert(qry<_eps)
    print('passed') 
# test linear scan
    print('Testing linear scan ... '.format(),end='')    
    step_size = 0.5
    origin = [0,-1,-1]
    axis = [1,-0.5, -1]
    steps = list(np.linspace(-3,3,num=7))
    scan = ScanGeometry.linear(steps, origin=origin, axis=axis, 
                               step_size=step_size)
    pos = scan.positions()
    # scan points
    r = asnumpy(pos['r'])
    # plot scan points in 3 projections
    lbl = [['y','z'],['z','x'],['x','y']]
    idx = [[1,2],[2,0],[0,1]]
    xlim = [-3,3]
    ylim = [-3,3]
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
    for i in range(3):
        ix, iy = idx[i]
        lblx, lbly = lbl[i]
        ax[i].set_xlim(xlim)
        ax[i].set_ylim(ylim)
        ax[i].plot(r[ix], r[iy], color='darkblue', marker='x', 
                   linestyle='-', markersize=5)
        ax[i].plot(r[ix,0], r[iy,0], color='darkblue', marker='o', 
                   linestyle='none', markersize=5)     
        ax[i].set_xlabel(lblx)
        ax[i].set_ylabel(lbly)
        ax[i].grid()
        ax[i].invert_xaxis()
    dstp = np.multiply(step_size*steps[0],axis)
    o = np.add(origin,dstp)
    title = 'origin: '+3*'{:g} '+'axis: '+3*'{:g} '+'step: {:g}'
    fig.suptitle(title.format(*o, *axis, step_size), fontsize=8)
    fig.tight_layout(w_pad=0.5)
    plt.show() 
    # assert value
    r0 = np.array([0,0,0])
    r_test = pos['tr'][4].r_to_loc(r0)
    rdiff = asnumpy(r_test) - np.array([0.5,-1.25, -1.5])
    qry = np.sum(rdiff**2)
    assert(qry<_eps)
    print('passed') 
    

if __name__ == "__main__":
    test_scan()




