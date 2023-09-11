# -*- coding: utf-8 -*-
"""
Defines basic surface elements with functions for event tracing.

Created on Mon Aug 14 18:48:33 2023
@author: Jan Saroun
"""


# TODO implement depth function as (i) projection along scandir, (ii) relative to
# a reference surface

# TODO: do it for both vi, vf vectors (path_in, path_out)

# TODO: benchmark cupy vs numpy
        
# TODO: write primitives also for: 
# plane, box, sphere, conus, ellipsoid, elliptic cylinder, hyper/parabola, toroid
# TODO: rewrite classes from the shapes module + add a general one
    
# TODO: develop visualization method (= plotContours from shapes)

import abc
# import numpy as np
import timeit
from cupyx.profiler import benchmark
import copy
import numpy as np

_gpu = False
def use_gpu(b):
    """Switch CuPy for GPU computing on/off."""
    global _gpu, cp
    if b:
        if _gpu:
            print('GPU mode is already switched on.')
        else:
            try:
                import cupy as cp
                _gpu = True
            except:
                print('Cannot import cupy.')        
                import numpy as cp
    else:
        if not _gpu:
            print('GPU mode is already switched off.')
        else:
            import numpy as cp
            _gpu = False


use_gpu(True)


_inf = 1e30

def _err_not_implemented():
    raise NotImplementedError('Subclass must implement abstract method.')
    

class Transform():
    """Coordinate transformation for surfaces and cells.
      
    Rotation is defined either by `angles`, or explicitely by rotation matrix.
    If `axes`=None, then `angles` correspond to the rotation axes x,y,z, 
    and the total rotation matrix is R = Rx.Ry.Rz. 
    Otherwise, `axes` defines the order of rotations. For example, 
    'axes' = [2,0,1] corresponds to R = Ry.Rx.Rz, and `angles` defines 
    corresponding rotation angles in the order [z,x,y].
    
    The translation is defined by the vectors `orig` and `rotctr`. 
    
    Transformation of a vector `r` from local to global coordinates is defined
    as:
    
    for order>0: r_glob = R.(r - rotctr) + orig + rotctr
    
    for order<0: r_glob = R.(r - orig - rotctr) + rotctr
    
    
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
        Order is translation + rotation (1) or rotation + translation (-1). 
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
            # inverse rotation
            self.imat = self.rmat.T
        # rotation is defined by axes and angles, Rot = R3.R2.R1
        else:
            A = axes['A']
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
        fmt = 'orig:' + 3*' {:g}'
        res = 'Transform('
        res += fmt.format(*self.orig[:,0])
        res += ', rotated: {}'.format(self.has_rotation)
        res += ')'
        return res
        
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
        return self._has_rotation

    @property
    def has_translation(self):
        return self._has_translation
    
    @property
    def is_identity(self):
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
        """Convert position to global coordinates."""
        if self.is_identity:
            return r
        elif not self._has_rotation:
            return r + self.orig
        else:
            return cp.dot(self.rmat,r) + self.orig
    
    def v_to_glob(self, v):
        """Convert direction vector to global coordinates."""
        if not self._has_rotation:
            return v
        else:
            return cp.dot(self.rmat,v)

    def r_to_loc(self, r):
        """Convert position to local coordinates."""
        if self.is_identity:
            return r
        elif not self._has_rotation:
            return r + self.iorig
        else:
            return cp.dot(self.imat,r) + self.iorig

    def v_to_loc(self, v):
        """Convert direction vector to local coordinates."""
        if not self._has_rotation:
            return v
        else:
            return cp.dot(self.imat,v)

class Surface():
    """Abstract base class for surface primitives."""
    
    def __init__(self, tr=None, inner=-1):
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
        
        """
        self._name = 'none' # id for debugging
        # local system: transformation w.r.t parent group
        if tr is None:
            self.trl = Transform() 
        else:
            self.trl = tr
        # global system: transformation w.r.t. root group.     
        self.trg = self.trl 
        self.inner = inner   # this is default by PHITS/MCNP convention
        self.nr = 0 # size of the event list
        self._reference = False # is it a reference surface?
    
    @property
    def reference(self):
        """Set True if this is a reference surface.
        
        A reference surface is used to define a coordinate system related 
        to the surface. Such a system has z-axis normal to the surface. The other
        axes depend on the surface shape, e.g. cylinder should have y||cylinder
        axis. The system is orthonormal with standard sign convention 
        (right-hand rule). 
        
        The surface coordinate system is used to calculate e.g. depth under 
        surface or components of a tensor field in other than cartesian 
        system. 
        
        """
        return self._reference
    
    @reference.setter
    def reference(self, value):
        self._reference = bool(value)
    
    def copy(self):
        """Deep copy if itself."""
        obj = type(self)()
        obj.__dict__ = copy.deepcopy(self.__dict__)
        return obj

    def connect(self, parent):
        """Connect coordinate system with the parent object."""
        self.trg = parent.trg.join(self.trg)
        #print('connected {}'.format(self))
        #print('       to {}'.format(parent))
    
    @abc.abstractmethod
    def set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.
        
        To be used for advance calculation of dependent arrays.
        
        Parameters
        ----------
        r : array_like
            Positions passed as (3,:) arrays.
        vi, vf : array_like
            Input and output velocities passed as (3,:) arrays.
        
        """
        _err_not_implemented()

    @abc.abstractmethod
    def inside(self, t=None, mask=None, side=-1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        t : ndarray (optional)
            Time shifts from the positions previously 
            defined by :meth:`~set_events`. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        path : str
            in|out for input|output directions.
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        """
        _err_not_implemented()
        
    @abc.abstractmethod
    def cross(self, path='in'):
        """Calculate cross-times with the surface.
        
        Given the position and velocity for each event, calculate:
        - Time to the entry and exit intersection with the surface
        
        The positions and velocities have to be previously defined
        by :meth:`set_events`
        
        Parameters
        ----------
        path : str
            in|out for input|output directions.
            
        Returns
        -------
        [xin, xout]
        
        xin, xout : list of input and output interesctions, respectively
            Each element is a pair of (time, mask), where `time` contains
            interestion times with the surface and `mask` indicates valid
            intersections.
        
        """
        _err_not_implemented()

    def evaluate(self):
        """Evaluate the object data to allow further calculations.
        
        Call this method before using `inside` and `cross` methods.
        
        """
        self.is_in = self.inside()
        
class Cylinder(Surface):     
    """Cylindric surface.
    
    Describes a cylinder along one of (x,y,z) axes.
    
    For axis=z, the equation is:
    (x − x0)^2 + (y − y0)^2 − R^2 = 0
    
    Parameters
    ----------
    R : float
        Radius.
    axis : str
        Axis direction, one of [x, y, z].       
    """
    
    _axes = {'x':[1,2], 'y':[0,2], 'z':[0,1]}
    def __init__(self, R=1.0, axis='z', **kwargs):         
        super().__init__(**kwargs)
        self.R = R
        self.axis = axis

    def __str__(self):
        res = 'Cylinder: R={:g}, trg={}'.format(self.R, self.trg)
        return res

    def set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed.

        see :func:`Surface.set_events`
        """
        self.nr = r.shape[1]
        r_loc = self.trg.r_to_loc(r)
        vi_loc = self.trg.v_to_loc(vi)
        vf_loc = self.trg.v_to_loc(vf)
        # velocity transfer
        self.q = vf_loc-vi_loc
        """Calculate what is needed to evaluate 'inside' mask and 
        intersections with the surface."""
        # local positions and velocities in x,z plane
        self.rc = r_loc[Cylinder._axes[self.axis],:]
        self.vic = vi_loc[Cylinder._axes[self.axis],:]
        self.vfc = vf_loc[Cylinder._axes[self.axis],:]
        # v^2 and r^2
        self.rsq = cp.sum(self.rc**2, axis=0)
        self.visq = cp.sum(self.vic**2, axis=0)
        self.vfsq = cp.sum(self.vfc**2, axis=0)
        # r.v
        self.rvi = cp.sum(self.rc*self.vic, axis=0)
        self.rvf = cp.sum(self.rc*self.vfc, axis=0)
        """ Calculate what is needed for a reference surface."""
        if self.reference:
            r0 = cp.sqrt(self.rsq)
            # surface coordinate system
            e0 = self.inner*r_loc[0]/r0
            e1 = cp.ones(self.nr)
            e2 = self.inner*r_loc[2]/r0
            nul = cp.zeros(self.nr)
            self.m_ref = cp.asarray([[e2, nul , -e0], 
                                    [nul, e1, nul], 
                                    [e0, nul, e2]])            
            # depth under the surface
            self._depth = self.inner*(r0 - self.R)
            # q in surface coordinates 
            self._qref = cp.einsum('ijk,jk->ik', self.m_ref, self.q)


    def inside(self, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        t : ndarray (optional)
            Time shifts from the positions previously 
            defined by :meth:`~set_events`. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        path : str
            in|out for input|output directions.
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        
        see :func:`Surface.inside`
        """    
        if t is None:
            rsq = self.rsq
        else:
            if path=='in':
                vc = self.vic
            else:
                vc = self.vfc
            rsq = cp.sum((self.rc+t*vc)**2, axis=0)
        D2 = rsq-self.R**2
        if mask is None:
            mask = cp.array(self.inner*side*D2 > 0, dtype=int)
        else:
            mask = mask*cp.array(self.inner*side*D2 > 0, dtype=int)
        return mask
    
    def cross(self, path='in'):
        """Calculate intersection times with the surface.
        
        see :meth:`Surface.cross`
        
        """
        if path=='in':
            rv = self.rvi
            vsq = self.visq
        else:
            rv = self.rvf
            vsq = self.vfsq
        D2 = rv**2 - vsq*(self.rsq-self.R**2)
        mask = cp.array(D2 > 0, dtype=int)
        D = cp.sqrt(cp.abs(D2))
        tin = (-rv - D)/vsq*mask + _inf*(1-mask)
        tout = (-rv + D)/vsq*mask + _inf*(1-mask)
        if self.inner<0:
            return [[tin], [tout]], [[mask], [mask]] 
        else:
            return [[tout], [tin]], [[mask], [mask]] 

    def depth(self):
        """Calculate depth under surface for each event.
        
        The meaning of "under surface" depends on the value of
        the `inner` attribute. For inner=-1, the depth is normal distance 
        calculated from the surface towards the cylinder axis.
        """
        if self._reference:
            result = self._depth
        else:
            return None

    def qref(self):
        """Get q-vectors in surface coordinates.
        
        q-vector is defined as velocity transfer, where velocities
        are defined as arguments when calling :meth:`set_events`.
        Defined only for reference surfaces, otherwise returns None.
        """
        if self._reference:
            result = self._qref
        else:
            return None

class Group(Surface):
# TODO rewrite for vi, vf   
# TODO: add functions for reference surface  
    def __init__(self, op='and', **kwargs):
        super().__init__(**kwargs)
        self.surfaces = []
        # initialize data for intersections
        self._isect = {}        
        for path in ['in', 'out']:
            self._isect[path] = {'times':[[],[]], 'mask':[[],[]]}       
        self._op = op.lower()
        if not self._op in ['and', 'or']:
            raise Exception('Invalid group operator: {}'.formatop())
    
    def __str__(self):
        res = 'Group: {}\n'.format(self.trg)
        for g in self.surfaces:
            res += '\t{}\n'.format(g)
        return res
    
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
        
        After adding all surfaces, call :meth:`set_events` and :meth:`evaluate`
        to make the object ready for use.
        
        Parameters
        ----------
        item : :class:`Surface`
            Instance of :class:`Surface` (can also be a Group) to be added.
        """    
        # create a copy of item and apply the chain of transformations
        item = item.copy()
        item.connect(self)
        self.surfaces.append(item)
    
    def set_events(self, r, vi, vf):
        """Set positions and velocities for all events to be processed."""
        self.nr = r.shape[1]
        for surf in self.surfaces:
            surf.set_events(r, vi, vf)

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
        for path in ['in', 'out']:
            """    
            The intersections are saved as a list of intersection times
            and corresponding mask which is 1 for intersections 
            to be **included**. Times are the flight times from the event 
            position to the intersection.  
            Assumes :func:`set_events` has been called to define the events.
            """
            
            xs_old = self._isect[path]['times'] # old intersections
            ms_old = self._isect[path]['mask'] # corresponding masks
            # for AND|OR:
            # keep old intersections only if inside|outside the new surface
            for i in range(2): # loop for input,output intersections 
                xs = xs_old[i]
                ms = ms_old[i]
                ns = len(xs)
                for j in range(ns):
                    # is old intersection inside|outside item?
                    q = item.inside(t=xs[j]+eps, mask=ms[j], side=sgn, path=path)
                    # update mask
                    ms[j] = ms[j]*q
            # for AND|OR:
            # add new intersections if inside|outside all of the other surfaces
            xs_new, ms_new = item.cross(path=path) # new intersections to be added
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
                        q1 = s.inside(t=xs[j]-eps, mask=ms[j], side=sgn, path=path)
                        q = q*q1
                    xs_old[i].append(xs[j])
                    ms_old[i].append(q)
        slist.append(item)

    def evaluate(self):     
        """Evaluate the group to allow further calculations.
        
        Call this method after adding all child surfaces or groups
        and after calling :meth:`set_events` on the top group.
        
        """
        iside = {-1:[0,1], 1:[1,0]}  # switch indices if self.inner=1
        # clear intersection data
        for path in ['in', 'out']:
            self._isect[path]['times'] = [[],[]]
            self._isect[path]['mask'] = [[],[]]
        added = [] # local variable for processed surfaces
        # collect intersections from all child surfaces
        for item in self.surfaces:
            # evaluate newly added surface or group
            item.evaluate()
#            if self.trg.orig[0]>0:
            #if self._name=='root':
#                print(self)
            self._add_xsections(item, added)
        # find events inside the group
        self.is_in = self.inside()
        
        # select only valid intersections:
        # - shift times for invalid intersections to _inf
        # - sort
        #if self.trg.orig[2]==2:
#        if self.trg.orig[0]>0:
#        if self._name=='root':
#            print(self)
        for path in ['in', 'out']:  # loop for incident|final directions
            times = cp.asarray(self._isect[path]['times'])
            mask = cp.asarray(self._isect[path]['mask'])
            for i in [0,1]: # loop for input|output interfaces
                t = times[i]*mask[i] + _inf*(1-mask[i])
                idx = cp.argsort(t,axis=0)
                # copy sorted arrays back to _isect
                # switch input/output if self.inner=1
                j = iside[self.inner][i]
                self._isect[path]['times'][j] = cp.take_along_axis(t, 
                                                                   idx, axis=0)
                self._isect[path]['mask'][j] = cp.take_along_axis(mask[i], 
                                                                  idx, axis=0)
        # calculate total path inside the group to the event position
        _tin, _tout = self._isect['in']['times']
        _min, _mout = self._isect['in']['mask']
        q_in = cp.array(_tin<0, dtype=int)*_min
        q_out = cp.array(_tout<0, dtype=int)*_mout
        self.path_in = cp.sum(_tout*q_out - _tin*q_in,axis=0)
        # calculate total path inside the group from the event position
        _tin, _tout = self._isect['out']['times']
        _min, _mout = self._isect['out']['mask']
        q_in = cp.array(_tin>0, dtype=int)*_min
        q_out = cp.array(_tout>0, dtype=int)*_mout
        self.path_out = cp.sum(_tout*q_out - _tin*q_in,axis=0)

# TODO remove
    def evaluate_old(self):     
        """Evaluate the group to allow further calculations.
        
        Call this method after adding all child surfaces or groups
        and after calling :meth:`set_events` on the top group.
        
        """
        if self.op=='or':
            sgn = -1 # OR operator
        else:
            sgn = 1 # AND operator 
        """
        Store intersections in `_tin` and `_tout` arrays.       
        The intersections are saved as a list of (time, mask) pairs, where mask 
        is 1 for intersections to be **included**, time is the flight time
        from the points pre-defined by :func:`set_events` to the intersection.
        """
        self._tin = []
        self._tout = []
        self._min = []
        self._mout = []
        added = [] # local variable for processed surfaces
        xs_old = [self._tin, self._tout] # old intersections
        ms_old = [self._min, self._mout] # corresponding masks
        # collect intersections from all child surfaces
        for item in self.surfaces:
            # evaluate newly added surface or group
            item.evaluate()
            # get intersections for the new surface
            # for AND|OR:
            # keep old intersections only if inside|outside the new surface
            for i in range(2):
                xs = xs_old[i]
                ms = ms_old[i]
                ns = len(xs)
                for j in range(ns):
                    x = xs[j]
                    # is old intersection (x) inside|outside item?
                    q = item.inside(t=x, mask=ms[j], side=sgn)
                    # update mask
                    ms[j] = ms[j]*q
            # for AND|OR:
            # add new intersections if inside|outside all of the other surfaces
            xs_new, ms_new = item.cross() # new intersections to be added
            for i in range(2):
                xs = xs_new[i]
                ms = ms_new[i]
                ns = len(xs)
                if not ns==len(ms):
                    raise Exception('Wrong lengths of intersection arrays')
                for j in range(ns):
                    x = xs[j]
                    q = ms[j]
                    for s in added:
                        # x is inside|outside s?
                        q1 = s.inside(t=xs[j], mask=ms[j], side=sgn)
                        q = q*q1
                    xs_old[i].append(xs[j])
                    ms_old[i].append(q) 
            added.append(item)
        # find events inside the group
        self.is_in = self.inside()
        
        # select only valid events
        m_in = cp.asarray(self._min)
        t = cp.asarray(self._tin)
        t_in = t*m_in + _inf*(1-m_in)
        idx = cp.argsort(t_in,axis=0)
        self._tin = cp.take_along_axis(t_in, idx, axis=0)
        self._min = cp.take_along_axis(m_in, idx, axis=0)
                
        m_out = cp.asarray(self._mout)
        t = cp.asarray(self._tout)
        t_out = t*m_out + _inf*(1-m_out)
        idx = cp.argsort(t_out,axis=0)
        self._tout = cp.take_along_axis(t_out, idx, axis=0)
        self._mout = cp.take_along_axis(m_out, idx, axis=0)
        
        # calculate path to the event position
        q_in = cp.array(self._tin<0, dtype=int)*self._min
        q_out = cp.array(self._tout<0, dtype=int)*self._mout
        self.path_in = cp.sum(self._tout*q_out - self._tin*q_in,axis=0)
        # calculate path from the event position
        q_in = cp.array(self._tin>0, dtype=int)*self._min
        q_out = cp.array(self._tout>0, dtype=int)*self._mout
        self.path_out = cp.sum(self._tout*q_out - self._tin*q_in,axis=0)
        # 
 
        
    def inside(self, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        t : ndarray (optional)
            Time shifts from the initial positions as 
            defined by :meth:`~set_events`. If None, use zero times.
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
        if self.op=='and':
            q = cp.ones(self.nr, dtype=int)
            for s in self.surfaces:
                if t is None:
                    q = q*s.is_in
                else:
                    is_in = s.inside(t=t, mask=mask, side=1, path=path)
                    q = q*is_in
        else:
            q = cp.zeros(self.nr, dtype=int)
            for s in self.surfaces:
                if t is None:
                    q = cp.bitwise_or(q,s.is_in)
                else:
                    is_in = s.inside(t=t, mask=mask, side=1, path=path)
                    q = cp.bitwise_or(q,is_in)
        # invert result if we check the outer side
        if side*self.inner>0:
            q = 1-q
        return q
   
    def cross(self, path='in'):
        """Calculate intersection times with the surface.
        
        see :meth:`Surface.cross`
        
        """ 
        return self._isect[path]['times'], self._isect[path]['mask']
    
    
#%% Tests
from matplotlib import pyplot as plt

class Bench():
    """Benchmark various operations with CuPy."""
    def slicing(self, v, slc=[2,0,1]):
        s = cp.random.permutation(slc)
        return v[s[0],:]
        #return v[:,s[0]]

    def sortfnc(self, a):
        idx = cp.argsort(a, axis=1)
        #b = cp.take_along_axis(a, idx, axis=1)
        return idx
        
    def bench_slicing(self):
        v = cp.random.random((3,1000000))
        #v = cp.random.random((100000,3))
        res =  benchmark(self.slicing, (v,), n_repeat=2000)
        return res 
    
    def bench_cross(self, gpu=True, n_repeat=200):
        n = 1000000
        obj = Cylinder()
        r = 2*cp.random.random((3,n))
        v = cp.random.random((3,n))
        v2 = v**2
        vn = cp.sqrt(cp.sum(v2, axis=0))
        v = v/vn
        time = None
        if gpu and _gpu:
            obj.set_events(r, v)
            res =  benchmark(obj.cross, (r,v), n_repeat=n_repeat)
            print(res)
            time = cp.average(res.gpu_times) + cp.average(res.cpu_times)
        else:
            glob = globals()
            glob.update({'r':r, 'v':v})
            res = timeit.timeit('obj.cross(r, v)', 
                          globals=glob, 
                          setup='obj = Cylinder(); obj.set_events(r, v)',
                          number=n_repeat)
            time = res/n_repeat
        return time
    
    def bench_sorting(self, gpu=True, n_repeat=200):
        n = 100000
        v = cp.random.random((n,5))
        time = None
        b = self.sortfnc(v)
        if gpu and _gpu:
            res =  benchmark(self.sortfnc, (v,), n_repeat=n_repeat)
            print(res)
            time = cp.average(res.gpu_times) + cp.average(res.cpu_times)
        else:
            glob = globals()
            glob.update({'v':v, 'self':self})
            res = timeit.timeit('self.sortfnc(v)', 
                          globals=glob, 
                          number=n_repeat)
            time = res/n_repeat
        return time, b

def test_bench_cross():
    bench = Bench()
    use_gpu(True)
    t1 = bench.bench_cross(gpu=True, n_repeat=200)
    use_gpu(False)
    t2 = bench.bench_cross(gpu=False, n_repeat=200)
    print('Bench cross:')
    print('GPU: {:g}'.format(t1))
    print('CPU: {:g}'.format(t2))
    print('gain: {:g}'.format(t2/t1))

def test_bench_sort():
    print('Bench sort:')
    bench = Bench()
    use_gpu(True)
    t1, b1 = bench.bench_sorting(gpu=True, n_repeat=200)
    use_gpu(False)
    t2, b2 = bench.bench_sorting(gpu=False, n_repeat=200)
    print('GPU: {:g}'.format(t1))
    print('CPU: {:g}'.format(t2))
    print('gain: {:g}'.format(t2/t1))
    print('b1:\n{}'.format(b1[0:5,:]))
    print('b2:\n{}'.format(b2[0:5,:]))
    
class Test():
    def __init__(self):
        self.root = None
        self.surf = None
        """
        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (6, 4),
                  'axes.labelsize': 'x-large',
                  'axes.titlesize':'x-large',
                  'xtick.labelsize':'large',
                  'ytick.labelsize':'large'}
        plt.rcParams.update(params)
        """
        theta = np.arange(-np.pi/2, np.pi/2, 0.01)        
        self.co =  np.cos(theta)
        self.si = np.sin(theta)

    
    def create_surface(self, R=1, inner=-1, **kwargs):
        """Create Surface.
        
        Assume Cylinder || y in (x,z) plane.
        
        Parameters
        ----------
        R : float
            Cylinder radius.
        **kwargs :
            Arguments used to create Transform object.
        
        """
        # assume Cylinder || y in (x,z) plane
        tr=Transform(**kwargs)        
        c = Cylinder(axis='y', R=R, inner=inner, tr=tr)
        return c

    def create_group(self, surfaces=[], groups=[], op='or', 
                     color=(0.5, 0.5, 0.5, 0.15), inner=-1, **kwargs):
        """Create Group with given list of surfaces ad sub-groups.
        
        Parameters
        ----------
        surfaces : list
            List of Arguments passed to  :meth:`create_surface`.
        groups : list
            List of Arguments passed to  :meth:`create_group`.  
        op : str
            Group operation, 'or'|'and'
        color : tuple
            Group color
        inner : int
            Define which side of the group surface is the inner one. 
            -1 is inside closed surfaces. It corresponds to the PHITS/MCNP 
            convention.            
        **kwargs :
            Arguments used to create Transform object.
        
        """
        #print('create group {}'.format(kwargs))
        root = Group(op=op, inner=inner, tr=Transform(**kwargs))
        # add color attribute
        root.color=color
        # first add all sub-groups
        for g in groups:
            grp = self.create_group(**g)
            root.add(grp)
        # than add surfaces
        for s in surfaces:
            c = self.create_surface(**s)
            root.add(c)
        #print('group id: {}'.format(id(root)))
        #print(root)
        return root

    def _add_group_limits(self, group, xlim, ylim):
        # assume (x,z) projection plane
        for el in group.surfaces:
            if isinstance(el, Group):
                xlim, ylim = self._add_group_limits(el, xlim, ylim)
            else: # assume cylinder
                o = el.trg.orig[0::2,0].get()
                xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
                xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
                ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
                ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        return xlim, ylim        
        
    def get_limits(self):
        xlim = [1e30, -1e30]
        ylim = [1e30, -1e30]
        self._add_group_limits(self.root, xlim, ylim)
        
        # make x,y ranges equal
        x0 = 0.5*(xlim[1] + xlim[0])
        y0 = 0.5*(ylim[1] + ylim[0])
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]
        scale = max(dx, dy)
        xlim[0] = x0 - 0.5*scale
        xlim[1] = x0 + 0.5*scale
        ylim[0] = y0 - 0.5*scale
        ylim[1] = y0 + 0.5*scale
  
        return xlim, ylim

    def set_events(self, n=10):
        xlim, ylim = self.get_limits() 
        x0 = 0.5*(xlim[1]+xlim[0])
        dx = xlim[1]-xlim[0]
        y0 = 0.5*(ylim[1]+ylim[0])
        dy = ylim[1]-ylim[0]
        
        rn = cp.random.random((3,n)) - 0.5
        x = dx*rn[0,:] + x0
        y = cp.zeros(n)
        z = dy*rn[2,:] + y0
        r = cp.array([x,y,z])
        #n=2
        #r = cp.array([[0, 0.45], [0, 0], [-2, -2]])
        vi = cp.array([cp.zeros(n), cp.zeros(n), cp.ones(n)])
        vf = cp.array([cp.ones(n), cp.zeros(n), cp.zeros(n)])
        self.r = r
        self.vi = vi
        self.vf = vf
        self.root.set_events(r, vi, vf)

    def set_events_test(self):
        #r = cp.array([[-2, 0, -1], [-3, 0, 0.5]]).T
        r = cp.array([[-1.5, 0, -0.5]]).T
        n = r.shape[1]
        vi = cp.array([cp.zeros(n), cp.zeros(n), cp.ones(n)])
        vf = cp.array([cp.ones(n), cp.zeros(n), cp.zeros(n)])
        self.r = r
        self.vi = vi
        self.vf = vf
        self.root.set_events(r, vi, vf)
        
    
    def plot_surface(self, ax, surf, color):  
        # assume Cylinder || y in (x,z) plane
        o = surf.trg.orig[0::2,0].get()
        x = surf.R*self.si + o[0]
        y = surf.R*self.co
        y2 = y + o[1]
        y1 = -y + o[1]
        ax.fill_between(x, y1, y2, facecolor=color)   
        
    def plot_group(self, ax, group, color):
        for el in group.surfaces:
            if isinstance(el, Group):
                cl = color
                if hasattr(el,'color'):
                    cl = el.color
                self.plot_group(ax, el, cl)
            else: 
                self. plot_surface(ax, el, color)
        
    def plot_cell(self, ax):
        gray = (0.5, 0.5, 0.5, 0.15)
        self.plot_group(ax, self.root, gray)     

    def plot_trace(self, ax):
        """Plot incident and output rays."""
        root = self.root
        lines = {'in':'b-','out':'r-'}
        sgn = {'in':1,'out':-1}
        # collect all incident rays
        for path in ['in','out']:
            rays = []
            _min,_mout = root._isect[path]['mask']
            _tin, _tout = root._isect[path]['times']
            isin = root.is_in
            if path=='in':
                q_in = cp.array(_tin<0, dtype=int)*_min
                q_out = cp.array(_tout<0, dtype=int)*_mout
            else:
                q_in = cp.array(_tin>0, dtype=int)*_min
                q_out = cp.array(_tout>0, dtype=int)*_mout
            for i in range(self.r.shape[1]):
                for j in range(len(_tin)):
                    if root.is_in[i]:
                        if path=='in':
                            if q_in[j,i]:
                                rin = (self.r[:,i] + self.vi[:,i]*_tin[j,i])[0::2]
                                if q_out[j,i]:
                                    rout = (self.r[:,i] + self.vi[:,i]*_tout[j,i])[0::2]
                                else:
                                    rout = self.r[:,i][0::2]
                                rays.append([rin,rout])
                        else:
                            if q_out[j,i]:
                                if q_in[j,i]:
                                    rin = (self.r[:,i] + self.vf[:,i]*_tin[j,i])[0::2]
                                else:
                                    rin = self.r[:,i][0::2]
                                rout = (self.r[:,i] + self.vf[:,i]*_tout[j,i])[0::2]
                                rays.append([rin,rout])
            if len(rays)>0:
                rays = cp.array(rays).get()
                cross_in = rays[:,0,:]
                cross_out = rays[:,1,:]
                ax.plot(cross_in[:,0], cross_in[:,1], 'b.', markersize=2)
                ax.plot(cross_out[:,0], cross_out[:,1], 'r.', markersize=2)
                for i in range(rays.shape[0]):
                    x = rays[i,:,0]
                    y = rays[i,:,1]
                    ax.plot(x, y, lines[path], linewidth=0.5)

    def plot(self, title='', trace=True):
        xlim, ylim = self.get_limits()
        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.plot_cell(ax)
        r0 = cp.asnumpy(self.r)
        if trace:
            m = cp.asnumpy(self.root.is_in)
        else:
            m = np.ones(r0.shape[1])
        r1 = r0[:,(m==0)]
        r2 = r0[:,(m>0)]
        ax.plot(r1[0], r1[2], color='gray', marker='x', linestyle='none', markersize=3)
        ax.plot(r2[0], r2[2], color='black', marker='x', linestyle='none', markersize=4)
        if trace:
            self.plot_trace(ax)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.grid()
        if title:
            plt.title(title, fontsize = 8)
        plt.show() 
    
    
    def group_moon(self):
        """Conjunction of inner and outer surfaces."""
        sf1 = {'R':1, 'inner':1, 'orig':[0.0, 0.0, 0.0]}           
        sf2 = {'R':1, 'orig':[0.5, 0.0, 0.0]}
        root = self.create_group(surfaces=[sf1, sf2], op='and')
        return root

    def group_hexagon(self, R1=1.5, R2=2,
                            op=['or','or','or','and', 'and'], 
                            inner=[1,1,1,-1, -1]):
        """Combine 6 surfaces in 3 groups + outer surface."""
        # hexagon parameter
        a = 2 
        co = np.cos(30*np.pi/180)
        si = 0.5
        # a pair of spheres
        sf1 = [{'R':R1, 'orig':[0.0, 0.0, 0.0]},
               {'R':R1, 'orig':[0.0, 0.0, 2*a*si]}
              ]
        
        gr = 3*[0]
        # hexagon from 3 pairs of spheres
        gr[0] = {'surfaces':sf1, 'orig':[a*co, 0., -a*si], 
               'op':op[0], 'inner':inner[0], 'color':(0.5, 0.0, 0.0, 0.15) }   
        
        gr[1] = {'surfaces':sf1, 'orig':[0, 0., a], 'angles':[0, -120, 0],
               'op':op[1], 'inner':inner[1], 'color':(0.0, 0.5, 0.0, 0.15)}

        gr[2] = {'surfaces':sf1, 'orig':[-a*co, 0., -a*si], 'angles':[0, -240, 0],
               'op':op[2], 'inner':inner[2], 'color':(0.0, 0.0, 0.5, 0.15)}
        
        gr2 = [{'groups':gr,'op':op[3], 'inner':inner[3]}]
        
        # sphere through hexagon corners
        sf2 = [{'R':R2, 'inner':inner[4], 'orig':[0.0, 0.0, 0.0]}]
        
        
        root = self.create_group(groups=gr2, surfaces=sf2, op=op[4])
        return root

    
    def test(self, plot=True, n=20):
        """Test tracing through groups of surfaces."""
        def get_title(op, inner):
            sgn = ['-','','']
            s = '{}[ '.format(sgn[inner[3]+1])
            for i in range(3):
                s += '{}(S{} {} S{})'.format(sgn[inner[i]+1], 2*i+1, op[i], 2*i+2)
                if i<2:
                    s += ' {} '.format(op[3])
            s += ' ] {} {}S7'.format(op[4], sgn[inner[4]+1])
            return s
        
        def test1(self):
            self.root = self.group_moon()
            self.set_events(n=n)
            #self.plot(trace=False)
            self.root.evaluate()
            if plot:
                self.plot()

        def test2(self):
            op=['or','or','or','and','and']; inner=[1,1,1,-1,-1]
            self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
            self.set_events(n=n)
            self.root.evaluate()
            if plot:
                self.plot(title=get_title(op, inner))

        def test3(self):
            op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, 1]
            self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
            self.set_events(n=n)
            self.root.evaluate()
            if plot:
                self.plot(title=get_title(op, inner))

        def test4(self):
            op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, -1]
            self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
            self.set_events(n=n)
            self.root.evaluate()
            if plot:
                self.plot(title=get_title(op, inner))

        def test5(self):
            op=['and','or','and', 'or','and']; inner=[-1,-1,-1,-1, -1]
            self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
            self.set_events(n=n)
            self.root.evaluate()
            if plot:
                self.plot(title=get_title(op, inner))

        def test0(self):
            #op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, -1]
            #op=['or','or','or','and','and']; inner=[1,1,1,-1,-1]
            op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, -1]
            self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)
            self.root._name='root'
            self.set_events_test()
            self.plot(trace=False)
            self.root.evaluate()
            if plot:
                self.plot(title=get_title(op, inner))
                
        
        my_tests = [test1, test2, test3, test4, test5]
        #my_tests = [test0]
        #my_tests = [test4, test5]
        for i in range(len(my_tests)):   
            print('Test {} '.format(i+1), end='')
            try:
                my_tests[i](self)
                print('passed')
            except Exception as e:
                print('failed')
                raise(e)
                #print(e)
        
        
#%% Run test
if __name__ == "__main__":
    test = Test()
    test.test(n=100, plot=True)





    

