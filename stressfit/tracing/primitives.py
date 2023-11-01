# -*- coding: utf-8 -*-
"""
Defines basic surface elements with functions for event tracing.

Created on Mon Aug 14 18:48:33 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""


# TODO implement depth function as (i) projection along scandir, (ii) relative to
# a reference surface

# TODO: benchmark cupy vs numpy
        
# TODO: write primitives also for: 
# plane, box, sphere, conus, ellipsoid, elliptic cylinder, hyper/parabola, toroid
# TODO: rewrite classes from the shapes module + add a general one

# TODO: make a separate module to handle cupy calls if cuda not available 
# (silent fallback to numpy + implement cp.asnumpy and cp.get)
    
import abc
# import numpy as np
import timeit
from cupyx.profiler import benchmark
import copy
import numpy as np
from .cuda import cp, asnumpy, has_gpu

_inf = 1e30

def _err_not_implemented(name):
    msg = 'Subclass must implement abstract method: {}.'.format(name)
    raise NotImplementedError(msg)
    

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
        _err_not_implemented('set_events')

    @abc.abstractmethod
    def inside(self, r=None, t=None, mask=None, side=-1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        r : ndarray (optional)
            Event positions in **global** coordinates. If None, use the events 
            defined by :meth:`~set_events`.
        t : ndarray (optional)
            Time shifts from the event positions. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        path : str
            in|out for input|output directions. Only used if t is defined.
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface.
        """
        _err_not_implemented('inside')
        
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
        _err_not_implemented('cross')

    @abc.abstractmethod
    def cal_depth(self):
        """Calculate depth under this surface for all events."""
        _err_not_implemented('cal_depth')

    @abc.abstractmethod    
    def cal_qref(self):
        """Calculate scattering vectors in surface coordinates.
        
        This is abstract method to be overriden by particular implementatin of 
        a primitive surface class.
        
        """
        _err_not_implemented('cal_qref')

    def evaluate(self):
        """Evaluate the object data to allow further calculations.
        
        Call this method before using `inside` and `cross` methods.
        
        """
        self.is_in = self.inside()


    def get_map(self, tr=None, depth=0.0, xlim=[-1,1], ylim=[-1,1], 
                npix=(500, 500)):
        """Create mask mapping the cut through the surface.
        
        By default, (x,z) cutting plane is assumed. Use Transform passed
        as the argument to cut through a different plane.
        
        Parameters
        ----------
        tr : Transform
            Optional root transformation which yields the cutting plane as 
            the (x,z) plane.
        depth : float
            y-coordinate of the cut.
        xlim : array_like(2)
            x-axis (abscisa) range.
        ylim : array_like(2)
            y-axis (dependent variable) range
        npix : tuple(2)
            number of pixels along x,y axes.
        
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
        m1 = self.inside(r=r1)
        m = asnumpy(m1.reshape(npix))
        return m
        
    
    
    
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

    def cal_depth(self):
        """Calculate depth under this surface for all events."""
        depth = None
        if self.nr>0:
            r0 = cp.sqrt(self.rsq)
            depth = self.inner*(r0 - self.R)
        return depth

    def cal_qref(self):
        """Calculate scattering vectors in surface coordinates.
        
        Surface coordinates are defined as
        
            - y || cylinder axis
            - z || radial coordinate of event position
        
        """
        qref = None
        if self.nr>0:
            # surface coordinate system
            r0 = cp.sqrt(self.rsq)
            e0 = self.inner*self.rc[0]/r0
            e1 = cp.ones(self.nr)
            e2 = self.inner*self.rc[1]/r0
            nul = cp.zeros(self.nr)
            m_ref = cp.asarray([[e2, nul , -e0], 
                                [nul, e1, nul], 
                                [e0, nul, e2]])
            qref = cp.einsum('ijk,jk->ik', m_ref, self.q)
        return qref

    def inside(self, r=None, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        r : ndarray (optional)
            Event positions in **global** coordinates. If None, use the events 
            defined by :meth:`~set_events`.
        t : ndarray (optional)
            Time shifts from the event positions. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        path : str
            in|out for input|output directions. Only used if t is defined.
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        
        see :func:`Surface.inside`
        """   
        if r is None:
            rc = self.rc
            rsq = self.rsq
        else:
            r_loc = self.trg.r_to_loc(r)
            rc = r_loc[Cylinder._axes[self.axis],:]
            rsq = cp.sum(rc**2, axis=0)
        if t is not None:
            if path=='in':
                vc = self.vic
            else:
                vc = self.vfc
            rsq = cp.sum((rc+t*vc)**2, axis=0)
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

class Group(Surface):  
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
        self.time_in = cp.sum(_tout*q_out - _tin*q_in,axis=0)
        # calculate total path inside the group from the event position
        _tin, _tout = self._isect['out']['times']
        _min, _mout = self._isect['out']['mask']
        q_in = cp.array(_tin>0, dtype=int)*_min
        q_out = cp.array(_tout>0, dtype=int)*_mout
        self.time_out = cp.sum(_tout*q_out - _tin*q_in,axis=0)
       
    def inside(self, r=None, t=None, mask=None, side=1, path='in'):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        r : ndarray (optional)
            Event positions in **global** coordinates. If None, use the events 
            defined by :meth:`~set_events`.
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
        if r is None:
            nr = self.nr
        else:
            nr = r.shape[1]
        if self.op=='and':
            q = cp.ones(nr, dtype=int)
            for s in self.surfaces:
                if (t is None) and (r is None):
                    q = q*s.is_in
                else:
                    is_in = s.inside(r=r, t=t, mask=mask, side=1, path=path)
                    q = q*is_in
        else:
            q = cp.zeros(nr, dtype=int)
            for s in self.surfaces:
                if (t is None) and (r is None):
                    q = cp.bitwise_or(q,s.is_in)
                else:
                    is_in = s.inside(r=r, t=t, mask=mask, side=1, path=path)
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
# TODO extract tests to separate module

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
        if gpu and has_gpu():
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
        if gpu and has_gpu():
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


