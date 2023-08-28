# -*- coding: utf-8 -*-
"""
Defines basic surface elements with functions for event tracing.

Created on Mon Aug 14 18:48:33 2023
@author: Jan Saroun
"""

# TODO: implement cross function for Group
# sort _xin, valid intersections to be at the begining
# sort _xout, valid intersections to be at the begining
# get valid index ranges
# TEST: we should get [in,out,in,out,in,out ...] sequence when sorted xin,xout 
# together.

# TODO: add complete transformations to primitives and groups

# TODO: implement path function

# implement depth function as (i) projection along scandir, (ii) relative to
# a reference surface

# TODO: do it for both vi, vf vectors (path_in, path_out)

# TODO: write tests for some simple groups and primitives

# TODO: write tests for a group of groups ... 

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
        self._has_rotation = abs(sp-3.0)<1e-10
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
        self._has_translation = s < 1e-10
            
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

# TODO pokr tady    # vyresit sekvenci pripojovani groups  
    
    def join(self, t):
        """Create new Transform which combines self and t."""
        if self.is_identity:
            obj = self
        else:            
            obj = Transform()
            obj.rmat = cp.dot(self.rmat, t.rmat)
            obj.imat = cp.dot(t.imat, self.imat)
            obj.orig = self.orig + cp.dot(self.rmat,t.orig)
            obj.iorig = t.iorig + cp.dot(t.imat,self.iorig)
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
    
    def __init__(self):
        """Initiate surface."""
        self.trl = Transform() # local system: transformation w.r.t parent group
        self.trg = self.trl # global system: transformation w.r.t. root group. 
        self.inner = -1   # this is default by PHITS/MCNP convention
        self.nr = 0 # size of the event list
    
    def copy(self):
        """Deep copy if itself."""
        obj = type(self)()
        obj.__dict__ = copy.deepcopy(self.__dict__)
        return obj

    def connect(self, parent):
        """Connect coordinate system with the parent object."""
        self.trg = parent.trg.join(self.trl)
    
    @abc.abstractmethod
    def set_events(self, r, v):
        """Set positions and velocities for all events to be processed.
        
        To be used for advance calculation of dependent arrays.
        
        Parameters
        ----------
        r, v : array_like
            Positions and velocities passed as (3,:) arrays.
        
        """
        _err_not_implemented()

    @abc.abstractmethod
    def inside(self, t=None, mask=None, side=-1):
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
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        """
        _err_not_implemented()
        
    @abc.abstractmethod
    def cross(self):
        """Calculate cross-times with the surface.
        
        Given the position and velocity for each event, calculate:
        - Time to the entry and exit intersection with the surface
        
        The positions and velocities have to be previously defined
        by :meth:`set_events`
            
        Returns
        -------
        [xin, xout]
        
        xin, xout : list of input and output interesctions, respectively
            Each element is a pair of (time, mask), where `time` contains
            interestion times with the surface and `mask` indicates valid
            intersections.
        
        """
        _err_not_implemented()


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
    def __init__(self, R=1.0, axis='z'):         
        super().__init__()
        self.R = R
        self.axis = axis

    def set_events(self, r, v):
        """Set positions and velocities for all events to be processed.

        see :func:`Surface.set_events`
        """
        self.nr = r.shape[1]
        r1 = self.trg.r_to_loc(r)
        v1 = self.trg.v_to_loc(v)
        self.rc = r1[Cylinder._axes[self.axis],:]
        self.vc = v1[Cylinder._axes[self.axis],:]
        self.rsq = cp.sum(self.rc**2, axis=0)
        self.vsq = cp.sum(self.vc**2, axis=0)
        self.rv = cp.sum(self.rc*self.vc, axis=0)

    def evaluate(self):
        """Evaluate the object data to allow further calculations.
        
        Call this method before using `inside` and `cross` methods.
        
        """
        self.is_in = self.inside()

    def inside(self, t=None, mask=None, side=1):
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
        
        Returns
        -------
        mask : array of int
            Mask (0,1) for valid events on given side of the surface. 
        
        see :func:`Surface.inside`
        """    
        if t is None:
            rsq = self.rsq
        else:
            rsq = cp.sum((self.rc+t*self.vc)**2, axis=0)
        D2 = rsq-self.R**2
        if mask is None:
            mask = cp.array(self.inner*side*D2 > 0, dtype=int)
        else:
            mask = mask*cp.array(self.inner*side*D2 > 0, dtype=int)
        return mask
    
    def cross(self):
        """Calculate intersection times with the surface.
        
        see :meth:`Surface.cross`
        
        """         
        D2 = self.rv**2 - self.vsq*(self.rsq-self.R**2)
        mask = cp.array(D2 > 0, dtype=int)
        D = cp.sqrt(cp.abs(D2))
        tin = (-self.rv - D)/self.vsq*mask + _inf*(1-mask)
        tout = (-self.rv + D)/self.vsq*mask + _inf*(1-mask)
        if self.inner<0:
            return [[tin], [tout]], [[mask], [mask]] 
        else:
            return [[tout], [tin]], [[mask], [mask]] 

class Group(Surface):
    def __init__(self, op='and'):
        super().__init__()
        self.surfaces = []
        self._tin = []   # input times
        self._tout = []  # mask for input times
        self._min = []   # output times
        self._mout = []  # mask for output times
        self._op = op.lower()
    
    @property
    def op(self):
        """Operator for grouping."""
        return self._op
            
    def connect(self, parent):
        """Connect coordinate system with the parent object."""
        super().connect(parent)
        for item in self.surfaces:
            item.connect(parent)            
    
    def add(self, item:Surface, inner=-1, tr=None):
        """Add a new Surface object to the group.
        
        After adding all surfaces, call :meth:`set_events` and :meth:`evaluate`
        to make the object ready for use.
        
        Parameters
        ----------
        item : :class:`Surface`
            Instance of :class:`Surface` (can also be a Group) to be added.
        inner : int
            (-1|1) indicates which side of the surface is the inner one. 
        tr : Transform
            Coordinate transformation relative to this group. See also 
            :class:`Transform`
        """
    
        # create a copy of item and apply the chain of transformations
        item = item.copy()
        item.inner = inner
        if tr is not None:
            item.trl = tr
        item.connect(self)
        self.surfaces.append(item)
    
    def set_events(self, r, v):
        """Set positions and velocities for all events to be processed."""
        self.nr = r.shape[1]
        for surf in self.surfaces:
            surf.set_events(r, v)

    def evaluate(self):
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
        
    def inside(self, t=None, mask=None, side=1):
        """Get mask for events on given surface side.
        
        Parameters
        ----------
        t : ndarray (optional)
            Time shifts from the initial positions as 
            defined by :meth:`~set_events`. If None, use zero times.
        mask : array of int (optional)
            Mask (0,1) for valid t elements.
        side : int
            Which side to check: 1:-1 stand for inner|outer sides.
            The inner side is defined by the `inner` attribute.
        
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
                    is_in = s.inside(t=t, mask=mask)
                    q = q*is_in
        else:
            q = cp.zeros(self.nr, dtype=int)
            for s in self.surfaces:
                if t is None:
                    q = cp.bitwise_or(q,s.is_in)
                else:
                    is_in = s.inside(t=t, mask=mask)
                    q = cp.bitwise_or(q,is_in)
        return q
   
    def cross(self):
        """Calculate intersection times with the surface.
        
        see :meth:`Surface.cross`
        
        """ 
        return [self._tin, self._tout], [self._min, self._mout] 
    

#%%
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
    
#%%

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

# test_bench_sort()

#%%
from matplotlib import pyplot as plt

class Test():
    def __init__(self):
        self.root = None
        self.surf = None
        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (6, 4),
                  'axes.labelsize': 'x-large',
                  'axes.titlesize':'x-large',
                  'xtick.labelsize':'large',
                  'ytick.labelsize':'large'}
        plt.rcParams.update(params)
    
    def create_group(self, R=[1., 1.5, 2.],
                     ctr=[[-1, 0., -1.],[0.5, 0., -0.5],[0.5, 0., 0.5]],
                     op='or'):
        root = Group(op=op)
        for i in range(len(R)):
            c = Cylinder(axis='y', R=R[i])
            root.add(c, tr=Transform(orig=ctr[i]))
        return root

    def define_cell(self, **kwargs):
        self.root = self.create_group(**kwargs)
        self.surf = self.root.surfaces
        
    def plot_cell(self,ax):
        theta = np.arange(-np.pi/2, np.pi/2, 0.01)
        gray = (0.5, 0.5, 0.5, 0.15)
        co =  np.cos(theta)
        si = np.sin(theta)
        for s in self.surf:
            o = s.trg.orig[0::2,0].get()
            x = s.R*si + o[0]
            y = s.R*co
            y2 = y + o[1]
            y1 = -y + o[1]
            ax.fill_between(x, y1, y2, facecolor=gray)

    def get_limits(self):
        xlim = [1e30, -1e30]
        ylim = [1e30, -1e30]
        for s in self.surf:
            o = s.trg.orig[0::2,0].get()
            xlim[0] = min(xlim[0],1.05*(-s.R+o[0]))
            xlim[1] = max(xlim[1],1.05*(s.R+o[0]))
            ylim[0] = min(ylim[0],1.05*(-s.R+o[1]))
            ylim[1] = max(ylim[1],1.05*(s.R+o[1]))
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
        v = cp.array([cp.ones(n), cp.zeros(n), cp.ones(n)])
        self.r = r
        self.v = v
        self.root.set_events(r, v)
        self.root.evaluate()
        

    def plot_trace(self, ax):
        rays = []
        root = self.root
        for i in range(self.r.shape[1]):
            for j in range(len(root._tin)):
                if root._min[j,i]:
                    rin = (self.r[:,i] + self.v[:,i]*root._tin[j,i])[0::2]
                    rout = (self.r[:,i] + self.v[:,i]*root._tout[j,i])[0::2]
                    rays.append([rin,rout])
        if len(rays)<=0:
            return
        rays = cp.array(rays).get()
        cross_in = rays[:,0,:]
        cross_out = rays[:,1,:]
        ax.plot(cross_in[:,0], cross_in[:,1], 'b.')
        ax.plot(cross_out[:,0], cross_out[:,1], 'r.')
        for i in range(rays.shape[0]):
            x = rays[i,:,0]
            y = rays[i,:,1]
            ax.plot(x, y, 'b-')


    def plot(self):
        xlim, ylim = self.get_limits()
        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.plot_cell(ax)
        r0 = cp.asnumpy(self.r)
        ax.plot(r0[0], r0[2], 'k+')
        self.plot_trace(ax)
        ax.grid()
        plt.show() 
    
    def test(self):
        par = {}
        par['R'] = [1., 1.3, 1.5, 0.5]
        par['ctr'] =[[-1, 0., -1.],
                     [0.5, 0., -0.5],
                     [0.5, 0., 0.5],
                     [-2.5, 0, -2.5]]
        
        self.define_cell(op='and', **par)
        self.set_events(n=20)
        self.plot()

#%%
        
test = Test()
test.test()



#%%
"""
tr =  Transform(orig=[5.0, 3., 1.], angles=[-45., 90., -30.])
rm = tr.rmat
o = tr.orig
r1 = cp.dot(rm, r) + o.reshape((3,1))
r1 = tr.v_to_glob(r)
r2 = tr.v_to_loc(r1)
print(r)
print(r1)
print(r2-r)
"""

