# -*- coding: utf-8 -*-
"""
Defines basic surface elements with functions for event tracing.

Created on Mon Aug 14 18:48:33 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

# TODO: unify set_events and evaluate, save memory ...
        
# TODO: write primitives also for: 
# plane, box, ellipsoid, elliptic cylinder, hyper/parabola, toroid
# TODO: rewrite classes from the shapes module + add a general one

    
import abc
import copy
import gc
import numpy as np
from scipy.spatial.transform import Rotation
from .cuda import cp, asnumpy

_inf = 1e30

# defined primitive surface shapes
_primitives = ['Cylinder','Plane', 'Sphere', 'Cone']

# all valid surface types including the predefined groups
_surfaces = _primitives + ['Box']

def is_primitive(obj):
    """Return true if obj is a primitive surface.
    
    Return false for composed (group) surfaces.
    """
    if isinstance(obj, Surface):
        s = obj.type_id
    else:
        s = obj
    return s in _primitives

def _err_not_implemented(name):
    msg = 'Subclass must implement abstract method: {}.'.format(name)
    raise NotImplementedError(msg)
    

def create_surface(shape='Cylinder', inner=-1, tr=None, **kwargs):
    """Create a Surface instance of given shape.
    
    Parameters
    ----------
    shape : str
        shape ID string
    inner : int
        Define which side of the surface is the inner one. 
        -1 is inside closed surfaces. It corresponds to the PHITS/MCNP 
        convention.
    tr : Transform
        Coordinate transformation relative to parent group. See also 
        :class:`Transform`
    **kwargs :
        Arguments specific to the shape passed to the constructor.
    
    """
    if not shape in _surfaces:
        raise Exception('Surface shape {} not defined'.format(shape))
    if shape=='Cylinder':
        c = Cylinder(inner=inner, tr=tr, **kwargs)
    elif shape=='Plane':
        c = Plane(inner=inner, tr=tr, **kwargs)
    elif shape=='Sphere':
        c = Sphere(inner=inner, tr=tr, **kwargs)
    elif shape=='Cone':
        c = Cone(inner=inner, tr=tr, **kwargs)
    elif shape=='Box':
        c = Box(inner=inner, tr=tr, **kwargs)
    return c

def create_transform(tr):
    if tr is None:
        tr = Transform() 
    elif isinstance(tr, dict):
        tr = Transform(**tr)
    elif isinstance(tr, Transform):
        pass
    else:
        print('tr = {}'.format(tr))
        raise Exception('Unknown coordinate transformation format.')
    return tr

class Transform():
    """Coordinate transformation for surfaces and cells.
      
    Rotation is defined either by `angles`, or explicitely by rotation matrix.
    If `axes`=None, then `angles` correspond to the rotation axes x,y,z, 
    and the total rotation matrix is R = Rx.Ry.Rz. 
    Otherwise, `axes` defines the order of rotations. For example, 
    'axes' = [2,0,1] corresponds to R = Ry.Rx.Rz, and `angles` defines 
    corresponding rotation angles in the order [z,x,y].
    
    The translation is defined by the vectors `orig` and `rotctr`. 
    
    Transformation of a vector `r` from local to global coordinates is:
    
    .. highlight:: python
    .. code-block:: python
    
        if order>0: 
            r_glob = R.(r - rotctr) + orig + rotctr
        else:
            r_glob = R.(r - orig - rotctr) + rotctr


    The inverse transformation (global to local) is:
    
    .. highlight:: python
    .. code-block:: python
     
         if order>0: 
             r_loc = R.T.(r - orig - rotctr) + rotctr 
         else:
             r_loc = R.T.(r - rotctr) + orig + rotctr

    Parameters
    ----------
    orig : array_like
        Translation vector.
    angles : array_like
        Rotation angles.
    rotctr : array_like
        Centre of rotation.
    axes : array_like, optional
        Axes indices.
    rmat : array_like, optional
        Rotation matrix.
    order : int
        The order of rotatio-translation transformations, see docstrings. 
    unit : str
        Angular units: deg or rad.
    
    """
    
    def __init__(self, orig=[0,0,0], angles=[0,0,0], rotctr=[0,0,0],
                 axes=None, rmat=None, order=1, unit='deg'):       
        # rotation is defined by angles corresponding to axes=[x,y,z],
        # Rot = Rx.Ry.Rz
        if axes is None:
            # rotation
            if rmat is not None:
                self.rmat = cp.asarray(rmat)
            elif angles is not None:
                self.rmat = self._cal_rmat(angles,unit=unit)
            else:
                self.rmat = cp.eye(3, dtype=float)
        # rotation is defined by axes and angles, Rot = R3.R2.R1
        else:
            A = axes
            # rotation
            self.rmat = self._cal_rmat(angles, unit=unit, axes=A[::-1])
        # inverse rotation
        self.imat = self.rmat.T
            
        # is there any rotation?
        sp = cp.sum(cp.diag(self.rmat))
        self._has_rotation = abs(float(sp)-3.0)>1e-10
        if not self._has_rotation:
            self.rmat = cp.eye(3, dtype=float)
            self.imat = cp.eye(3, dtype=float)
    
    
        # evaluate translation
        if orig is None:
            orig = [0.,0.,0.]
        if rotctr is None:
            rotctr = [0.,0.,0.]
        orig = cp.asarray(orig).reshape((3,1))
        rotctr = cp.asarray(rotctr).reshape((3,1))
        if order>0:
            # translation
            self.orig = orig + rotctr - cp.dot(self.rmat,rotctr)
            # inverse translation
            self.iorig = rotctr - cp.dot(self.imat, orig+rotctr)
        else:
            # translation
            self.orig = rotctr - cp.dot(self.rmat, orig+rotctr)
            # inverse translation
            self.iorig = orig + rotctr - cp.dot(self.imat,rotctr)
        
        # is there any translation?
        s = cp.sum(cp.abs(self.orig)) + cp.sum(cp.abs(self.iorig))
        self._has_translation = float(s) > 1e-10
    
    def __str__(self):
        res = 'Transform('
        if self.has_translation:
            fmt = 'orig: [' + 3*' {:g}'+']'
            o = cp.round(self.orig[:,0],10)
            res += fmt.format(*o)
        if self.has_rotation:
            if self.has_translation:
                res += ', '
            angles = self._get_angles()
            fmt = 'angles: [' + 3*' {:g}'+']'
            o = np.round(angles,6)
            res += fmt.format(*o)
        res += ')'
        return res
    
    def _get_angles(self):
        """Calculate xyz Euler angles from global rotation matrix."""
        r = asnumpy(self.rmat)
        r =  Rotation.from_matrix(r)
        angles = r.as_euler("xyz",degrees=True)
        return angles
    
    def _cal_rmat(self, angles, unit='deg', axes=[0,1,2]):
        deg = cp.pi/180
        idx = [[1,2],[2,0],[0,1]]
        if unit=='deg':
            u = deg
        else:
            u = 1
        R = 3*[None]
        for i in range(3):
            R[i] = cp.eye(3, dtype=float)
            c = cp.cos(angles[i]*u)
            s = cp.sin(angles[i]*u)
            ix, iy = idx[axes[i]]
            R[i][ix,ix] = c
            R[i][iy,iy] = c
            R[i][ix,iy] = -s
            R[i][iy,ix] = s
        rmat = cp.dot(R[0],cp.dot(R[1],R[2]))
        return rmat
    
    @property
    def has_rotation(self):
        """Rotation is included."""
        return self._has_rotation

    @property
    def has_translation(self):
        """Translation is included."""
        return self._has_translation
    
    @property
    def is_identity(self):
        """Identity transformation."""
        return not (self._has_rotation or self._has_translation)

    def join(self, t):
        """Create new Transform which combines self and t."""
        if self.is_identity:
            obj = t
        else:            
            obj = Transform()
            obj.rmat = cp.dot(self.rmat, t.rmat)
            obj.imat = cp.dot(t.imat, self.imat)
            obj.orig = self.orig + cp.dot(self.rmat,t.orig)
            obj.iorig = t.iorig + cp.dot(t.imat,self.iorig)
            # is there any translation?
            s = cp.sum(cp.abs(obj.orig)) + cp.sum(cp.abs(obj.iorig))
            obj._has_translation = float(s) > 1e-10
            sp = cp.sum(cp.diag(obj.rmat))
            obj._has_rotation = abs(float(sp)-3.0)>1e-10
        return obj
    
    def r_to_glob(self, r):
        """Transform position to global coordinates."""
        # ensure correct handling of single r-vectors
        rr = r
        if len(r.shape)<2: rr = r.reshape(3,-1)
        if self.is_identity:
            res = cp.asarray(rr)
        elif not self._has_rotation:
            res = cp.asarray(rr) + self.orig
        else:
            res = cp.dot(self.rmat,rr) + self.orig
        return res.reshape(r.shape)
    
    def v_to_glob(self, v):
        """Transform direction vector to global coordinates."""
        if not self._has_rotation:
            return cp.asarray(v)
        else:
            return cp.dot(self.rmat,cp.asarray(v))

    def r_to_loc(self, r):
        """Transform position to local coordinates."""
        # ensure correct handling of single r-vectors
        rr = r
        if len(r.shape)<2: rr = r.reshape(3,-1)
        if self.is_identity:
            res = cp.asarray(r)
        elif not self._has_rotation:
            res = cp.asarray(rr) + self.iorig
        else:
            res = cp.dot(self.imat,cp.asarray(rr)) + self.iorig
        return res.reshape(r.shape)

    def v_to_loc(self, v):
        """Transform direction vector to local coordinates."""
        if not self._has_rotation:
            return cp.asarray(v)
        else:
            return cp.dot(self.imat,cp.asarray(v))

#%% Surface primitives

class Surface():
    """Abstract base class for surface primitives."""
    
    def __init__(self, tr=None, inner=-1, level=0, is_reference=False):
        """Initiate surface.
        
        Parameters
        ----------
        tr : Transform
            Coordinate transformation relative to parent group. See also 
            :class:`Transform`
        inner : int
            Define which side of the surface is the inner one. 
            -1 is inside closed surfaces. It corresponds to the PHITS/MCNP 
            convention.
        level : int
            Level in the group hierarchy.
        is_reference : bool
            Set true to make it a reference, which allows to calculate
            principal direction and depth coordinates relative to the surface.
            Only primitive surfaces can be defined as a reference.
        
        """
        # accuracy limit [mm] for boundary determination
        # keep at zero to avoid other problems !!
        self._eps = 0.0 
        # storage for temporary data
        self._tmp = {}
        # storage for r,vi,vf in global coordinates
        self._glob = []
        # transformation relative to the parent 
        self.trl = create_transform(tr)
        # transformation relative to the global reference frame
        self.trg = self.trl 
        self.inner = inner    # denotes the inner side of the surface 
        self.is_reference = is_reference   # set as reference? 
        self._level = level   # level in the group hierarchy.
        
        
# TODO: following is needed only for a reference surface        
        # Depth under the surface for each event - only reference surfaces
        self._depth = None 
        # Q-vectors in principal reference frame of the surface
        self._qref = None
        # Q-vectors in local coordinates
        self.q = None

    @property
    def is_reference(self):
        """Is a reference surface."""
        return self._reference
    
    @is_reference.setter
    def is_reference(self, value:bool):
        if self.type_id in _primitives:
            self._reference = value
        else:
            self._reference = False

    @property
    def type_id(self):
        """Return class name."""
        return type(self).__name__
    
    @property
    def depth(self):
        """Return depth under this surface for all events."""
        res = None
        if self._reference:
            res = self._depth
        return res

    @property
    def qref(self):
        """Return q-vectors in principal reference frame of the surface.
        
        Only available for a reference surface.
        """
        res = None
        if self._reference:
            res = self._qref
        return res

    def get_events(self, path):
        """Return event data in local coordinates.
        
        Use data previously defined by :meth:`evaluate`.
        
        Parameters
        ----------
        path : str
            in|out for incident and exit rays
            
        Returns
        -------
        r, v : array_like
            Positions and velocities for given path. If path is not defined, 
            then v is a list with both incident and exit velocities.
        
        """
        idx = {'in':1, 'out':2}
        if path in idx:
            iv = idx[path]
        else:
            iv = [1,2]
        if 'events' in self._tmp:
            r = self._tmp['events'][0]
            v = self._tmp['events'][iv]
        else:
            # check that there are valid data in _glob
            q1 = len(self._glob)==3
            q2 = [isinstance(a,cp.ndarray) for a in self._glob]
            if q1 and all(q2):
                ev = [0,0,0]
                ev[0] = self.trg.r_to_loc(self._glob[0])
                for i in [1,2]:
                    ev[i] = self.trg.v_to_loc(self._glob[i])
                r = ev[0]
                v = ev[iv]
            else:
                raise Exception('No event data stored.')
        if isinstance(v,list):
            return [r] + v
        else:
            return r,v

    def _param_str(self):
        """Return class-specific parameters as a string."""
        return ''
    
    def __str__(self):
        """Return string representation of the object."""
        res = '{}: '.format(self.type_id)
        p = self._param_str()
        if p:
            res += ' {}'.format(p)
        if not self.trg.is_identity:
            if p:
                res += ','
            res += ' {}'.format(self.trg)
        return res    
    
    def copy(self):
        """Deep copy if itself."""
        obj = type(self)()
        obj.__dict__ = copy.deepcopy(self.__dict__)
        return obj

    def connect(self, parent):
        """Connect coordinate system with the parent object."""
        self.trg = parent.trg.join(self.trg)
    
    def cleanup(self):
        """Clean temporary data storage to save memory."""
        self._tmp.clear()

    def cal_qref(self, U, q):
        """Calculate scattering vectors in given coordinates.
        
        Parameters
        ----------
        U : ndarray
            Unitary matrix with basis vectors given as rows. 
            For planes, a single matrix (3,3) is provided. Otherwise,
            the reference frame depends on each event position
            and the matrix has the shape (3,3,:), where the last index 
            denotes respective events. 
        
        q: ndarray
            A (3,:) array of q-vectors in local coordinates.
        
        """
        rank = len(U.shape)
        if rank == 2:
            qref = cp.dot(U, q)
        elif rank==3:
            qref = cp.einsum('ijk,jk->ik',U, q)
        else:
            qref = None
        return qref
        
    def get_map(self, tr=None, depth=0.0, xlim=[-1,1], ylim=[-1,1], 
                npix=(500, 500)):
        """Create mask mapping the cut through the surface.
        
        By default, (x,z) cutting plane is assumed. Use Transform passed
        as the argument to cut through a different plane.
        
        Parameters
        ----------
        tr : Transform
            Optional root transformation which transforms x,y map
            coordinates to the cut [x,y,z] coordinates: 
            [x, y, z]_cut = tr.r_to_loc([x, depth, y])
        depth : float
            y-coordinate of the cut.
        xlim : array_like(2)
            x-axis (abscisa) range.
        ylim : array_like(2)
            y-axis (dependent variable) range
        npix : tuple(2)
            number of pixels along x,y axes.
        
        Returns
        -------
        array
            Pixel map as NumPy array (converted from CuPy if applicable).
        
        """        
        if tr is not None:
            t = tr
        else:
            t = Transform()
        # define coordinates of map bins
        dx = -np.subtract(*xlim)/npix[0]
        dy = -np.subtract(*ylim)/npix[1]
        x = cp.linspace(*xlim, num=npix[0], endpoint=False, dtype=float)+0.5*dx
        y = cp.linspace(*ylim, num=npix[1], endpoint=False, dtype=float)+0.5*dy          
        xx = cp.resize(x,npix)
        yy = cp.resize(y,(npix[1],npix[0])).T
        zz = depth*cp.ones(npix)
        r = cp.asarray([xx.reshape(-1),zz.reshape(-1),yy.reshape(-1)])
        r1 = t.r_to_loc(r)
        m1 = self.is_inside(r=r1)
        m = asnumpy(m1.reshape(npix))
        return m

    def is_inside(self, r=None, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        If `r` or `t` is not specified, use the pre-calculated data 
        for events specified when calling :meth:`evaluate`.
        
        Parameters
        ----------
        r : ndarray (optional)
            Event positions in **global** coordinates. If None, use the events 
            defined by :meth:`~evaluate`.
        t : ndarray (optional)
            Time shifts from the event positions, using provided velocities. 
            If t=None or v=None, no shift is applied.
        mask : array of int (optional)
            Additional mask (0,1) to be applied on the result.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute of the surface.
        path : str
            in|out for input|output directions. Only used if t is defined.

        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface.
        """   
        
        rloc = None
        if (r is not None):
            rloc = self.trg.r_to_loc(r)
        if (t is not None):
            ev = self.get_events(path)
            if rloc is None:
                rloc = ev[0]
            rloc = rloc + ev[1]*t
        # NOTE: it is OK to pass rloc=None to _inside
        isin = self._inside(rloc=rloc)
        if mask is not None:
            isin = isin*mask
        if side>0:
            return isin
        else:
            return 1-isin

    def _cross_times(self):
        """Get intersection times and masks for this surface.
        
        Uses local events pre-calculated by the :meth:`evaluate` method.
        
        Returns
        -------        
        The dict structure X[path][var][dir] = [xin, xout], where
        
        path : str
            'in' or 'out' for incident and output ray directions, respectively
        var : str
            'time' for cross-times and 'mask' for (0|1) mask marking the valid 
            cross-sections
        dir : int
            0|1 for inward|outward cross points
        
        xin, xout : array_like
            The array(s) of crossing times/masks (depending on the `var` value). 
            A list of multiple arrays is possible for each xin, xout if the 
            surface allows for multiple cross points.
        """
        # calculate intersections
        res = {}        
        for path in ['in', 'out']:
            res[path] = {}
            tx, mx = self._cross(path)
            res[path]['time'] = tx
            res[path]['mask'] = mx
        return res

    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        To be used for advance calculation of temporary arrays.
        By default, it converts event data to local reference frame and saves 
        results in the _tmp['events'] field. 
        
        Parameters
        ----------
        r : array_like
            Positions passed as (3,:) arrays.
        vi, vf : array_like
            Input and output velocities passed as (3,:) arrays.
        
        """
        rloc = self.trg.r_to_loc(r)
        viloc = self.trg.v_to_loc(vi)
        vfloc = self.trg.v_to_loc(vf)
        self._tmp['events'] = [rloc, viloc, vfloc]
        self._glob = [r, vi, vf]

    
    def evaluate(self, r, vi, vf):
        """Evaluate the surface for given event coordinates.
        
        Call this method before using `is_inside` and `cross_times` methods.
        
        Parameters
        ----------
        r : array_like
            Positions passed as (3,:) arrays.
        vi, vf : array_like
            Input and output velocities passed as (3,:) arrays.
        
        """
        self.cleanup()
        self._set_events(r, vi, vf)
        #isin = self._inside()
        #cross = self._cross_times()
        if self.is_reference:
            self._cal_reference()


# Abstract methods to be implemented by child classes        
    @abc.abstractmethod
    def _cal_reference(self):
        """Calculate data required for a reference surface.
        
        This method must define following fields for each event, which 
        has been defined in :meth:`_set_events`:
        
        - `_depth`: Depth under the reference surface.
        
        - `_qref`: Q-vectors in principal surface coordinates. 
          
        The surface basis vectors are expressed as rows of the matrix U
        in local coordinates. The matrix U is used to calculate principal 
        strain/stress components for given surface symmetry. 
        
        The definition of the principal surface coordinates depends 
        on the surface shape. Usually the Z-component
        should be along the nearest surface normal, Y denotes
        axial and X circumferential (hoop) components. 

        Parameters
        ----------
        rloc : ndarray (optional)
               Event positions in **local** coordinates.

        """
        _err_not_implemented('_cal_reference')

    @abc.abstractmethod
    def _inside(self, rloc=None):
        """Get mask for events on given surface side.
        
        Abstract method to be implemented by each primitive surface class.
        
        Parameters
        ----------
        rloc : ndarray (optional)
               Event positions in **local** coordinates. If not defined,
               the method should use values pre-calculated in 
               :meth:`_set_events`
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface.
        """
        _err_not_implemented('_inside')

    @abc.abstractmethod
    def _cross(self, path):
        """Calculate cross-times with the surface.

        Abstract method to be implemented by each primitive surface class.
        
        Use :meth:`get_events(path)` to get the positions and velocities
        in **local** coordinates:
        
        r, v = self.get_events(path)
        
        For each event (columns of r,v), calculate:
            
        - Time to the entry and exit intersection with the surface
        - Mask invalid intersections (no intersection, corresponding times
          are set to +- inf)
        
        Parameters
        ----------
        path : str
            in|out for incident and exit rays
            
        Returns
        -------
        [tin, tout], [min, mout]
        
        tin, tout : list
            Lists of inward and outward interesction times as arrays.
            Both tin, tout can include multiple arrays if there are more 
            than one input/output cross-section.
        min, mout : list
            Corresponding masks indicating valid intersections.
        """
        _err_not_implemented('_cross')

class Plane(Surface):     
    """A Plane.
    
    Describes a cylinder along one of (x,y,z) axes.
    
    The plane is defined by equation:
    p[0]*x + p[1]*y + p[2]*z − d = 0
    
    
    The surface-related basis vectors (X,Y,Z) are defined as follows:

    - Z = plane normal
    - X = y x Z if y<>Z, else X = x
    - Y = Z x X 
        
    where (x,y,z) is the local basis.    
    
    Parameters
    ----------
    p : array_like
        Normal vector to the plane.
    d : float
        Distance from origin.    
    """
    
    def __init__(self, p=[0,1,0], d=0.0, **kwargs):         
        super().__init__(**kwargs)
        n = cp.array(p)
        a = cp.linalg.norm(n)
        # plane normal and distance, normalized
        self.n = n/a
        self.d = float(d/a)

    #override
    def _param_str(self):
        fmt = 'p=['+2*'{:g},'+'{:g}] d={:g}'
        res = fmt.format(*self.n, self.d)
        return res

    #override
    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        See docstrings of :meth:`Surface._set_events`.
        """
        super()._set_events(r, vi, vf)
        rloc = self._tmp['events'][0]
        self._tmp['rn'] = cp.dot(self.n, rloc) - self.d

    #override
    def _cal_reference(self):
        # calculate depth under reference surface
        viloc, vfloc = self._tmp['events'][1:]
        rn = self._tmp['rn']
        self._depth = self.inner*rn
        eps = 1e-8  
        y = cp.array([0,1,0], dtype=float)
        x = cp.array([1,0,0], dtype=float)
        Z = self.n
        X = cp.cross(y,Z)
        if cp.linalg.norm(Z)<eps:
            X = x
        Y = cp.cross(Z, X)
        U = cp.asarray([X,Y,Z])
        q = vfloc - viloc
        # calculate q in reference surface coordinates
        self._qref = self.cal_qref(U, q)
        
    #override
    def _inside(self, rloc=None):  
        if rloc is None:
            rn = self._tmp['rn']
        else:
            rn = cp.dot(self.n, rloc) - self.d    
        mask = cp.array(self.inner*rn > self._eps, dtype=int)
        return mask
    
    #override
    def _cross(self, path):
        _eps = 1e-20
        v = self.get_events(path)[1] 
        rn = self._tmp['rn']
        vn = cp.dot(self.n, v)
        t = -rn/(vn+_eps)
        # mask for valid cross-sections
        m = cp.array(abs(vn) > _eps, dtype=int)        
        # mask for crossing inwards        
        m_in = m*cp.array(self.inner*vn > 0, dtype=int)
        m_out = m*(1 - m_in)
        # times for crossing inwards
        tin =  _inf*m_out + t*m_in
        # times for crossing outwards
        tout = _inf*m_in + t*m_out
        return [[tin], [tout]], [[m_in], [m_out]] 
    

class Sphere(Surface):     
    """A Sphere.
    
    Describes a spherical surface.
    
    The surface basis (X,Y,Z) is expressed in terms of the 
    local coordinates (x,y,z). 
    
    The surface-related basis vectors (X,Y,Z) are defined as follows:

    - Z = surface normal
    - X = y x Z if y<>Z, else X = x
    - Y = Z x X
        
    where (x,y,z) is the local basis.        
    
    Parameters
    ----------
    R : float
        Radius.
    **kwargs :
        Other arguments passed to the parent Surface constructor. 

    """
    
    def __init__(self, R=1.0, **kwargs):         
        super().__init__(**kwargs)
        self.R = abs(R)

    def _param_str(self):
        res = 'R={:g}'.format(self.R)
        return res

    def _C(self, r=None):
        if r is None:
            C = self._tmp['C']
        else:
            rsq = cp.sum(r**2, axis=0)
            C = rsq - self.R**2
        return C

    #override
    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        See docstrings of :meth:`Surface._set_events`.
        """
        super()._set_events(r, vi, vf)
        rloc, viloc, vfloc = self._tmp['events']
        rsq = cp.sum(rloc**2, axis=0)
        self._tmp['rsq'] = rsq
        self._tmp['C'] = rsq - self.R**2
        self._tmp['rv']  = {'in':cp.sum(rloc*viloc, axis=0), 
                            'out':cp.sum(rloc*vfloc, axis=0)}
        self._tmp['vsq'] = {'in':cp.sum(viloc**2, axis=0), 
                            'out':cp.sum(vfloc**2, axis=0)}

    #override
    def _cal_reference(self):
        eps = 1e-10
        rsq = self._tmp['rsq']
        rloc, viloc, vfloc = self._tmp['events']
        r0 = cp.sqrt(rsq)
        self._depth = self.inner*(r0 - self.R)
        m = cp.array(r0>eps, dtype=int)
        x = cp.array([1,0,0],dtype=float)
        y = cp.array([0,1,0],dtype=float)
        Z = self.inner*rloc/(r0+(1-m)*eps)
        X = m*cp.cross(y,Z) + (1-m)*x
        Y = cp.cross(Z,X)
        U = cp.asarray([X,Y,Z])
        q = vfloc - viloc
        # calculate q in reference surface coordinates
        self._qref = self.cal_qref(U, q)

    #override
    def  _inside(self, rloc=None): 
        C = self._C(r=rloc)
        mask = cp.array(self.inner*C > self._eps**2, dtype=int)
        return mask
    
    #override
    def _cross(self, path):
        r, v = self.get_events(path)
        C = self._tmp['C']
        rv = self._tmp['rv'][path]
        vsq = self._tmp['vsq'][path]
        D2 = rv**2 - vsq*C
        mask = cp.array(D2 > 0, dtype=int)
        D = cp.sqrt(cp.abs(D2))
        tin = (-rv - D)/vsq*mask + _inf*(1-mask)
        tout = (-rv + D)/vsq*mask + _inf*(1-mask)
        if self.inner<0:
            return [[tin], [tout]], [[mask], [mask]] 
        else:
            return [[tout], [tin]], [[mask], [mask]] 
    
class Cylinder(Sphere):     
    """Cylindric surface.
    
    Describes a cylinder along one of (x,y,z) axes.
    
    For axis=z, the equation is:
    x**2 + y**2 − R**2 = 0
    
    The surface-related basis vectors (X,Y,Z) are defined as follows:
    
    - Y = cylinder axis
    - Z = surface normal along (x,y,0), or (1,0,0) if not defined
    - X = Y x Z if y<>Z, else X = x 
    
    where (x,y,z) is the local basis.    
    
    Parameters
    ----------
    R : float
        Radius.
    axis : str
        Axis direction, one of [x, y, z]. 
    **kwargs :
        Other parameters passed to the parent Sphere constructor.
    """
    
    
    def __init__(self, R=1.0, axis='z', **kwargs):
        _rot = {'x':[0,90,0], 'y':[90,0,0], 'z':[0,0,0]}
        ax = axis.lower()
        if not ax in _rot:
            raise Exception('Invalid axis: {}'.format(ax))
        tr = Transform(angles=_rot[ax])
        if 'tr' in kwargs:
            tr0 = create_transform(kwargs['tr'])
            tr = tr0.join(tr)
            del kwargs['tr']
        super().__init__(R=R, tr=tr, **kwargs)
        self.axis = axis.lower()
            

    def _param_str(self):
        res = 'R={:g}, axis={}'.format(self.R, self.axis)
        return res

    #override
    def _C(self, r=None):
        if r is None:
            C = self._tmp['C']
        else:
            rc = r[0:2,:]
            rsq = cp.sum(rc**2, axis=0)
            C = rsq - self.R**2
        return C

    #override
    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        See docstrings of :meth:`Surface._set_events`.
        """
        super()._set_events(r, vi, vf)
        rloc, viloc, vfloc = self._tmp['events']
        vin = viloc[0:2,:]
        vout = vfloc[0:2,:]
        rc = rloc[0:2,:]
        rsq = cp.sum(rc**2, axis=0)
        self._tmp['rsq'] = rsq
        self._tmp['C'] = rsq - self.R**2
        self._tmp['rv']  = {'in':cp.sum(rc*vin, axis=0), 
                            'out':cp.sum(rc*vout, axis=0)}
        self._tmp['vsq'] = {'in':cp.sum(vin**2, axis=0), 
                            'out':cp.sum(vout**2, axis=0)}

    #override       
    def _cal_reference(self):
        """Set the reference frame related to this surface.
        
        See docstring of :meth:`Surface.set_reference`
        """
        eps = 1e-10
        rloc, viloc, vfloc = self._tmp['events']
        nr = rloc.shape[1]
        one = cp.ones(nr)
        nul = cp.zeros(nr)
        rsq = self._tmp['rsq']
        r0 = cp.sqrt(rsq)
        self._depth = self.inner*(r0 - self.R)   
        x = cp.asarray([one,nul],dtype=float)
        m = cp.array(r0>eps, dtype=int)
        ez = self.inner*rloc[0:2,:]/(r0+(1-m)*eps)
        Z = m*ez + (1-m)*x
        U = cp.asarray([[-Z[1], Z[0], nul], 
                        [nul ,  nul , one],
                        [Z[0],  Z[1], nul]
                       ])
        q = vfloc - viloc
        # calculate q in reference surface coordinates
        self._qref = self.cal_qref(U, q)

class Cone(Cylinder):     
    """A Plane.
    
    Describes a cone along one of (x,y,z) axes.
    
    The cone equation is
    
    sqrt(x**2 + y**2) = abs(R*(z − a))
    
    for a cone axis='z'. Other directions and positions can be defined
    through **kwargs as the Transform parameter `tr` of the :class:`Surface`
    constructor. 
    
    The surface-related basis vectors (X,Y,Z) are defined as follows:

    - Z = nearest surface normal
    - Y = along nearest surface axis
    - X = circular (hoop) component
    
    Parameters
    ----------
    a : float
        Position of the appex on the cone axis. 
    R: float
       Cone radius at a unit distance from apex. Equals to the tangent of 
       the cone half-angle.
    axis : str
        Cone axis, one of x, y or z.
    **kwargs :
        Other parameters passed to the parent Surface constructor. 
    """
    
    def __init__(self, a=0.0, R=1, axis='z', **kwargs):         
        super().__init__(R=R, axis=axis, **kwargs)
        self.a = a

    def _param_str(self):
        fmt = 'a={:g}, R/{}={:g}, axis={}'
        res = fmt.format(self.a, self.axis, self.R, self.axis)
        return res
    
    #override
    def _C(self, r=None):
        if r is None:
            C = self._tmp['C']
        else:
            C = super()._C(r=r)
            ra = r[2,:] - self.a 
            C = C + self.R**2*(1-ra**2)
        return C
    
    #override
    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        See docstrings of :meth:`Surface._set_events`.
        """
        super()._set_events(r, vi, vf)
        rloc, viloc, vfloc = self._tmp['events']
        rsq = self._tmp['rsq']
        ra = rloc[2,:] - self.a 
        va = {'in':viloc[2,:],'out':vfloc[2,:]}
        # transform rsq etc. from cylinder to cone shape
        # then we can re-use the cylinder equations
        self._tmp['C'] = rsq - self.R**2*ra**2
        for path in ['in','out']:
            rv = self._tmp['rv'][path]
            vsq = self._tmp['vsq'][path]
            self._tmp['rv'][path]  = rv - self.R**2*ra*va[path]
            self._tmp['vsq'][path] = vsq - self.R**2*va[path]**2 

    #override       
    def _cal_reference(self):
        """Set the reference frame related to this surface.
        
        See docstring of :meth:`Surface.set_reference`
        """
        eps = 1e-10
        rloc, viloc, vfloc = self._tmp['events']
        nr = rloc.shape[1]
        ra = rloc[2,:] - self.a
        # cone angles cos and sin
        co = 1/np.sqrt(1+self.R**2)
        si = self.R*co
        # radial component
        x = cp.linalg.norm(rloc[0:2,:], axis=0)
        # axial component
        z = cp.abs(ra)
        self._depth = self.inner*(x*co - z*si) # OK, tested
        
        one = cp.ones(nr,dtype=float)
        nul = cp.zeros(nr,dtype=float)
        # mask for out-of-axis points
        mx = cp.array(x>eps, dtype=int)
        # azimuthal angle of the radial component
        # avoid division-by-zero errors
        # we could use arctan2, but this should be fatser ...
        if cp.sum(mx)<nr:
            cosia = mx*rloc[0:2,:]/(x+(1-mx)*eps) + (1-mx)*cp.array([one, nul])
        else:
            cosia = rloc[0:2,:]/x
        # sign of z-coordinate
        zsn = cp.sign(rloc[2,:])
        # principal axes for in-plane coordinates (r[1]=0)
          # z = cp.array([co*one,nul,-si*zsn]) # normal to surface
          # y = cp.array([si*one,nul,co*zsn]) # along surface axis
          # x = cp.array([nul,one,nul]) # hoop 
        # rotate by the azimuthal angle around cone axis
        Z = [co*cosia[0,:], -co*cosia[1,:], -si*zsn] # normal to surface
        Y = [si*cosia[0,:], -si*cosia[1,:], co*zsn]  # along surface axis
        X = [cosia[1,:], cosia[0,:], nul ]           # hoop 
        U = cp.asarray([X, Y, Z])
        q = vfloc - viloc
        # calculate q in reference surface coordinates
        self._qref = self.cal_qref(U, q)

#%% Group of surfaces

class Group(Surface):
    """Group of surface primitives or other surface groups.
    
    Parameters
    ----------
    op : str
        Operator for grouping, [and|or]
    surfaces : list
        List of Surface objects.
    groups : list
        List of Group objects.
    color : tuple
        Color specification
    **kwargs : 
        Other arguments passed to the parent Surface constructor. 
        See docstring for :class:`Surface`. 
    
    """
    
    def __init__(self, op='and', surfaces=[], groups=[], 
                 color=(0.5, 0.5, 0.5, 0.15), **kwargs):
        super().__init__(**kwargs)
        self._op = op.lower()
        if not self._op in ['and', 'or']:
            raise Exception('Invalid group operator: {}'.format(op))
        self.color = color
        # initialize data for intersections
        self._isect = {}        
        for path in ['in', 'out']:
            self._isect[path] = {'time':[[],[]], 'mask':[[],[]]}       
        self.surfaces = []
        # first add all sub-groups
        # TODO : put everything in the surfaces argument, no need to distinguish ... 
        for g in groups:
            if isinstance(g, Group):
                grp = g.copy()
            else:
                grp = Group(level=self._level+1,**g)
            #print('Create {}'.format(grp))
            self.add(grp)
        # then add surfaces
        for s in surfaces:
            if isinstance(s, Surface):
                c = s.copy()
            else:
                c = create_surface(level=self._level+1,**s)
            self.add(c)
        self.is_in = None

    def _isect_clean(self):
        for path in ['in', 'out']:
            self._isect[path]['time'] = [[],[]]
            self._isect[path]['mask'] = [[],[]]

    def cleanup(self):
        """Clean temporary data storage to save memory."""
        # clean temporary data in all child surfaces
        for s in self.surfaces:
            s.cleanup()
        if self._level>0:
            # clear intersection data except for the root
            self._isect_clean()
        else:
            # for root group: clean references to r,vi,vk in all childs
            for surf in self.surfaces:
                surf._glob = [0,0,0]
            # call garbage collector by the root group
            gc.collect()

    def __str__(self):
        fmt = self._level*'\t'+'{}\n'
        res = fmt.format(super().__str__())
        for g in self.surfaces:
            if isinstance(g,Group):
                res += g.__str__()
            else:
                fmt = g._level*'\t'+'{}\n'
                res += fmt.format(g)
        return res
    
    def is_member(self, surface:Surface):
        """Check if given Surface object is member of this surface group."""
        ans = (surface in self.surfaces)
        if not ans:
            for s in self.surfaces:
                if isinstance(s, Group):
                    ans = s.is_member(surface)
                else:
                    ans = (surface == s)
                if ans:
                    break       
    @property
    def op(self):
        """Operator for grouping."""
        return self._op
            
    def connect(self, parent):
        """Connect coordinate system with the parent object."""
        super().connect(parent)
        for item in self.surfaces:
            item.connect(parent)            
    
    def add(self, item:Surface):
        """Add a new Surface object to the group.
        
        After adding all surfaces, call :meth:`evaluate`
        to make the object ready for use.
        
        Parameters
        ----------
        item : :class:`Surface`
            Instance of :class:`Surface` (can also be a Group) to be added.
        """    
        # create a copy of item and apply the chain of transformations
        item.connect(self)
#        if isinstance(item, Group):
#            print('Add {} to {}\n{}'.format(item.type_id, self.type_id, item))
        self.surfaces.append(item)
    
    def _set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed."""
        self.nr = r.shape[1]
        for surf in self.surfaces:
            surf._set_events(r, vi, vf)

# TODO
    def _add_xsections(self, item, slist):
        """Add intersections for a new surface object.
        
        Parameters
        ----------
        item : Surface
            New Surface object to search for interfaces.
        slist : list
            List of Surfaces allready processed
            
        For a new Surface object (item):
            
        - calculate intersections with surfaces in the given list
        - combine interfaces (remove intersections which disapeared due to 
          adding operation - AND/OR)
        - add new intersections in the _isect field
        - add item to the list of already processed surfaces (slist)
        """
        eps = 1e-10 # small shift to avoid problems at intersections of surfaces
        signs = {'or':-1, 'and':1}
        sgn = signs[self.op] 
        # loop for incident|final ray directions
        
        cross = item._cross_times()
        
        """
        if isinstance(item, Group):
            rr = item.surfaces[0]._glob[0]
            ct = cross['in']['time']
            cm = cross['in']['mask']
            print(item)
            fmt = 3*'{:g} '+', '+2*'{:g} '+', '+2*'{:g} '
            print(fmt.format(*rr[:,0],*ct[0][:].T[0],*cm[0][:].T[0]))
            print('')
        """     
        
#        ERROR
        # TODO error - does not remember cross times for a group
        for path in ['in', 'out']:
            """    
            The intersections are saved as a list of intersection times
            and corresponding mask which is 1 for intersections 
            to be **included**. Times are the flight times from the event 
            position to the intersection.  
            """        
            xs_old = self._isect[path]['time'] # old intersections
            ms_old = self._isect[path]['mask'] # corresponding masks
            # for AND|OR:
            # keep old intersections only if inside|outside the new surface
            for i in range(2): # loop for input,output intersections 
                xs = xs_old[i]
                ms = ms_old[i]
                ns = len(xs)
                for j in range(ns):
                    # is old intersection inside|outside item?
                    # update mask
                    q = item.is_inside(t=xs[j]+eps, mask=ms[j], side=sgn, 
                                           path=path)
                    ms[j] = ms[j]*q
            # for AND|OR:
            # add new intersections if inside|outside all of the other surfaces
            xs_new, ms_new = cross[path].values() # new intersections to be added
            for i in range(2):
                xs = xs_new[i]
                ms = ms_new[i]
                ns = len(xs)
                if not ns==len(ms):
                    raise Exception('Wrong lengths of intersection arrays')
                for j in range(ns):
                    q = ms[j]
                    for s in slist:
                        # x is inside|outside s?
                        q1 = s.is_inside(t=xs[j]-eps, mask=ms[j], side=sgn, 
                                         path=path)
                        q = q*q1
                    xs_old[i].append(xs[j])
                    ms_old[i].append(q)
        slist.append(item)

    def evaluate(self, r, vi, vf):     
        """Evaluate the group to allow further calculations.
        
        Call this method after adding all child surfaces or groups
        and after calling :meth:`_set_events` on the top group.
        
        """
        iside = {-1:[0,1], 1:[1,0]}  # switch indices if self.inner=1
        self.cleanup()
        self._set_events(r, vi, vf)
        # clear intersection data
        self._isect_clean()
        added = [] # local variable for processed surfaces
        # collect intersections from all child surfaces
        if self.op=='and':
            isin = cp.ones(self.nr, dtype=int)
        else:
            isin = cp.zeros(self.nr, dtype=int)
        for item in self.surfaces:
          # evaluate newly added surface or group
          # NOTE: this will call evaluate + _add_xsections on all item's childs 
            #print(item)
            item.evaluate(r, vi, vf)
            mask = item.is_inside()
            if self.op=='and':
                isin = isin*mask
            else:
                isin = cp.array(isin+mask>0,dtype=int)
            self._add_xsections(item, added)
            item.cleanup() # cleanup to save memory
        self._tmp['inside'] = isin
        
        # select only valid intersections:
        # - shift times for invalid intersections to _inf
        # - sort
        #if self.trg.orig[2]==2:
#        if self.trg.orig[0]>0:
#        if self._name=='root':
#            print(self)
        
        ms = {'in':[0,0], 'out':[0,0]}
        xs = {'in':[0,0], 'out':[0,0]}
        for path in ['in', 'out']:  # loop for incident|final directions
            times = cp.asarray(self._isect[path]['time'])
            mask = cp.asarray(self._isect[path]['mask'])
            
            
            for i in [0,1]: # loop for input|output interfaces
                t = times[i]*mask[i] + _inf*(1-mask[i])
                idx = cp.argsort(t,axis=0)
                # copy sorted arrays back to _isect
                # switch input/output if self.inner=1
                j = iside[self.inner][i]
                xs[path][j] = list(cp.take_along_axis(t,idx, axis=0))
                ms[path][j] = list(cp.take_along_axis(mask[i], idx, axis=0))
                #self._isect[path]['time'][j] = cp.take_along_axis(t,idx, axis=0)
                #self._isect[path]['mask'][j] = cp.take_along_axis(mask[i], idx, axis=0)
        
        # accept only rays with some intersections
        #"""
        self._isect_clean() 
        for path in ['in', 'out']:
            nx = len(ms[path][0])
            for j in range(nx):
                q1 = ms[path][0][j].any()
                q2 = ms[path][1][j].any()
                if q1 and q2:
                #if True:
                    for i in [0,1]:
                        self._isect[path]['mask'][i].append(ms[path][i][j])
                        self._isect[path]['time'][i].append(xs[path][i][j])
        #"""    

        
# TODO create RootGroup(Group) and move following there  ?
        # calculate total path inside the group to the event position
        if self._level == 0:
            # calculate total path inwards
            _tin, _tout = cp.asarray(self._isect['in']['time'])
            _min, _mout = cp.asarray(self._isect['in']['mask'])
            q_in = cp.array(_tin<0, dtype=int)*_min
            q_out = cp.array(_tout<0, dtype=int)*_mout
            self.time_in = cp.sum(_tout*q_out - _tin*q_in,axis=0)
            # calculate total path outwards
            _tin, _tout = cp.asarray(self._isect['out']['time'])
            _min, _mout = cp.asarray(self._isect['out']['mask'])
            q_in = cp.array(_tin>0, dtype=int)*_min
            q_out = cp.array(_tout>0, dtype=int)*_mout
            self.time_out = cp.sum(_tout*q_out - _tin*q_in,axis=0)
            self.is_in = isin
               
    def _inside(self): 
        if 'inside' in self._tmp:
            # use _tmp
            q = self._tmp['inside']
        else:
            q = self.is_in
        return q      

       
    def is_inside(self, r=None, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        r : ndarray (optional)
            Event positions in **global** coordinates. If None, use the events 
            defined by :meth:`~evaluate`.
        t : ndarray (optional)
            Time shifts from the initial positions as 
            defined by :meth:`~evaluate`. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1|-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        path : str
            in|out for input|output directions.
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        """       
        if (t is None) and (r is None):
           q = self._inside()
        else:
            if r is None:
                nr = self.nr
            else:
                nr = r.shape[1]
            # this allows to call is_inside for t<>0 after cleaning surfaces
            if self.op=='and':
                q = cp.ones(nr, dtype=int)
                for s in self.surfaces:
                    q = q*s.is_inside(r=r, t=t, mask=mask, side=1, path=path)
            else:
                q = cp.zeros(nr, dtype=int)
                for s in self.surfaces:
                    is_in = s.is_inside(r=r, t=t, mask=mask, side=1, path=path)
                    q = cp.bitwise_or(q,is_in)
        # invert result if we check the outer side
        if side*self.inner>0:
            q = 1-q
        return q

    
    def _cross(self, path):
        """Calculate intersection times with the surface.
        
        see :meth:`Surface._cross`
        
        """ 
        return self._isect[path]['time'], self._isect[path]['mask']
    

class Box(Group):     
    """A rectangual box, composed of planes.  
    
    Parameters
    ----------
    size : array_like
        Box dimensions along x,y,z    
    """
    
    def __init__(self, size=[1,1,1], **kwargs):  
        self.size = size
        px1 = {'shape':'Plane', 'inner':1, 'p':[1,0,0],'d':-0.5*size[0]} 
        px2 = {'shape':'Plane', 'inner':-1, 'p':[1,0,0],'d':0.5*size[0]}
        py1 = {'shape':'Plane', 'inner':1, 'p':[0,1,0],'d':-0.5*size[1]} 
        py2 = {'shape':'Plane', 'inner':-1, 'p':[0,1,0],'d':0.5*size[1]}
        pz1 = {'shape':'Plane', 'inner':1, 'p':[0,0,1],'d':-0.5*size[2]} 
        pz2 = {'shape':'Plane', 'inner':-1, 'p':[0,0,1],'d':0.5*size[2]}
        gx = {'op':'and','surfaces':[px1, px2]}
        gy = {'op':'and','surfaces':[py1, py2]}
        gz = {'op':'and','surfaces':[pz1, pz2]}
        super().__init__(op='and', groups=[gx,gy,gz], **kwargs)
        
    def _param_str(self):
        fmt = 'size=['+2*'{:g},'+'{:g}]'
        res = fmt.format(*self.size)
        return res
