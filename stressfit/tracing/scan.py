# -*- coding: utf-8 -*-
"""
Scan geometry module.

Created on Sat Sep 16 15:36:02 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

from .cuda import cp
from .primitives import Transform

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
    
    Parameters
    ----------
    shape : int
        Scan shape, use ScanGeometry.LINEAR or CIRCULAR.
    steps : array_like
        Step positions in [mm] or [deg], depending on scan shape.
    origin : array_like
        Coordinates of the scan origin (step=0)
    axis : array_like
        Scan direction or rotation axis, depending on scan shape.
    step_size : float
        Step size factor (multiplies the values of `steps`) 
    origin_rot : array_like
        Origin of the rotation axis (a point on the rotation axis)    
    tr : Transform
        Initial sample positioning transformation.
    
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
                 step_size=1.0, origin_rot=[0,0,0], tr=None):
        if not shape in [0,1]:
            raise Exception('Unsupported scan shape: {}'.format(shape))
        self._shape = shape
        self._steps = cp.array(steps, dtype=float)
        self._origin = cp.array(origin, dtype=float)
        self._step_size = step_size
        a = cp.array(axis, dtype=float)
        self._axis = a/cp.linalg.norm(a) # normalize direction axis
        self._origin_rot = cp.array(origin_rot, dtype=float)
        if isinstance(tr, Transform):
            self._tr = tr
        else:
            self._tr = Transform()

    @classmethod
    def linear(cls, steps, origin=[0,0,0], axis=[0,1,0], step_size=1.0, tr=None):
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
        tr : Transform
            Initial sample positioning transformation.
        
        """
        c = cls(ScanGeometry.LINEAR, steps, origin=origin, axis=axis, 
                step_size=step_size, tr=tr)
        return c

    @classmethod
    def circular(cls, steps, origin=[1,0,0], axis=[0,1,0], origin_rot=[0,0,0],
                 step_size=1.0, tr=None):
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
        tr : Transform
            Initial sample positioning transformation.
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

    @property
    def transform(self):
        """Initial sample positioning transformation."""
        return self._tr
    
    @transform.setter
    def transform(self, value):
        if isinstance(value, Transform):
            self._tr = value
        elif isinstance(value, dict):
            self._tr = Transform(**value)

    @property
    def scan_orig(self, local=True):
        """Scan origin coordinates.
        
        Parameters
        ----------
        local : bool
            If true, return postion in cell root coordinates. 
            Otherwise, use the `transform` property (initial sample position)
            for conversion to global (lab) coordinates.
        """
        if local:
            return self._origin
        else:
            return self._tr.r_to_glob(self._origin)

    @property
    def rot_orig(self, local=True):
        """Rotation axis origin coordinates.
        
        Parameters
        ----------
        local : bool
            If true, express postion in cell root coordinates. 
            Otherwise, use the `transform` property (initial sample position)
            for conversion to global (lab) coordinates.
        """
        if local:
            return self._origin_rot
        else:
            return self._tr.r_to_glob(self._origin_rot)

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
            tr = Transform(rmat=R, orig=-self._origin, rotctr=self._origin_rot, 
                           order=1)
            #if not self._tr.is_identity:
               # tr = self._tr.join(tr)
            #    tr = tr.join(self._tr)
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
            #if not self._tr.is_identity:
            #    tr = self._tr.join(tr)
            trans.append(tr)
        return {'r': r.T, 'tr':trans}

    def positions(self):
        """Calculate transformations which generate scan positions.
        
        Returns
        -------
        dict
            r : array_like
                Scan positions in global (lab) coordinates, shape = (3,:).
            tr : list
                A list of Transform objects, which generate scan positions and 
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
        if self._shape == ScanGeometry.LINEAR:
            res = self._positions_lin()
        elif self._shape == ScanGeometry.CIRCULAR:
            res = self._positions_rot()
        else:
            res = None
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




