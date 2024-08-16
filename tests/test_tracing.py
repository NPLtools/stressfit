# -*- coding: utf-8 -*-
"""
Tests of stressfit.tracing package.

Created on Mon Sep 25 15:33:39 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

#%% Test objects
import abc
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from stressfit.tracing import options
options['gpu'] = True
from stressfit.tracing.cuda import cp, asnumpy, has_gpu, to_numpy
from stressfit.tracing.cells import Cell, Extinction, Material
from stressfit.tracing.primitives import Transform, Cylinder, Plane, Sphere, ECylinder
from stressfit.tracing.primitives import Group, RootGroup
from stressfit.tracing.events import StrainSampling
from stressfit.tracing.scan import ScanGeometry
from stressfit.tracing.integrals import GaugeIntegral
from stressfit.dataio import load_data
import pandas as pd

_run_tests = ['scan','primitives', 'integrals']
_run_tests = ['scan']
_run_tests = ['bench']
_run_tests = ['integrals']
_run_tests = ['primitives']
#_run_tests = []

# TODO test PolygonBar 
# TODO test with variable velocities
    
def test_transform():
    print('\ntest_transform ... ',end='')
    deg = np.pi/180
    L1 = 2
    L2 = 0.5
    dx1 = 1
    dx2 = 2
    ang1 = 30
    ang2 = 15
    t2 = Transform(orig=[L1,0,0], angles=[0,ang1,0], rotctr=[dx1,0,0])
    t1 = Transform(orig=[L2,0,0], angles=[0,ang2,0], rotctr=[dx2,0,0]) 
    r0 = cp.zeros(3)
    #r2 = t2.r_to_loc(r0)
    #r1 = t1.r_to_loc(r2)   
    tt = t2.join(t1)
    rr = tt.r_to_loc(r0)  
    #fmt = '{}: '+3*'{:g} '+', {:g}'
    #print(fmt.format('r0',*r0, cp.linalg.norm(r0)))
    #print(fmt.format('r2',*r2, cp.linalg.norm(r2)))
    #print(fmt.format('r1',*r1, cp.linalg.norm(r1)))
    #print(fmt.format('rr',*rr, cp.linalg.norm(rr)))
    x1 = (L1+dx1)*cp.cos(ang1*deg)-dx1
    y1 = (L1+dx1)*cp.sin(ang1*deg)
    arot = ang2*deg + np.arctan2(y1,x1+L2+dx2)
    R = np.sqrt((x1+L2+dx2)**2+y1**2)
    x2 = R*cp.cos(arot)-dx2
    y2 = R*cp.sin(arot)
    #print(fmt.format('cal1',-x1,0,-y1,np.sqrt(x1**2+y1**2)))
    #print(fmt.format('cal2',-x2,0,-y2,np.sqrt(x2**2+y2**2)))
    diff = cp.linalg.norm(cp.array([-x2-rr[0], -y2-rr[2]]))
    #print(diff)
    assert(diff<1e-10)
    print('OK')
    
def test_scan():
    """Test circular and linear scan settings."""
    test_transform()
    print('\ntest_scan:')
    _eps = 1e-15
# test circular scan
    step_size = 1
    steps = cp.linspace(-30,120,num=11)
    origin = [2,0,2]
    axis = [0,1,0]
    origin_rot = [-1,0,0]
    print('Circular scan ... ',end='')
    scan = ScanGeometry.circular(steps,origin=origin, axis=axis, 
                                 origin_rot=origin_rot, step_size=step_size)
    scan.transform = {'orig':[0,0,0], 'angles':[0, -45, 0], 'order':1}
    # test points
    pos = scan.positions()
    r = asnumpy(pos['r'])
    # ki, kf vectors, assume ki || z, kf || x
    x_ax = np.array([1,0,0])
    z_ax = np.array([0,0,1])
    a = {'ki':[], 'kf':[]}
    for i in range(scan.nsteps):
        x_axt = asnumpy(pos['tr'][i].v_to_loc(x_ax))
        z_axt = asnumpy(pos['tr'][i].v_to_loc(z_ax))
        kf = np.array([r[:,i], r[:,i]+x_axt]).T
        ki = np.array([r[:,i]-z_axt, r[:,i]]).T
        a['kf'].append(kf)
        a['ki'].append(ki)
    # plot scan points with attached ki, kf vectors
    xlim = [-4,4]
    ylim = [-4,4]
    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.plot(r[0], r[2], color='black', marker='x', linestyle='none', markersize=3)    
    ax.plot(origin[0], origin[2], color='red', marker='o', linestyle='none', 
            markersize=5)
    ax.plot(origin_rot[0], origin_rot[2], color='red', marker='^', 
            linestyle='none', markersize=5)
    for i in range(scan.nsteps):
        ax.plot(a['kf'][i][0], a['kf'][i][2], color='red', linestyle='-', 
                linewidth=0.5)
        ax.plot(a['ki'][i][0], a['ki'][i][2], color='blue', linestyle='-', 
                linewidth=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.grid()
    title = 'origin: '+3*'{:g} '
    title += 'axis: '+3*'{:g} '
    title += 'rot_origin: '+3*'{:g} '
    title += 'step: {:g}'
    plt.title(title.format(*origin, *axis, *origin_rot, step_size), fontsize=8)    
    #ax.invert_xaxis()
    plt.show()     
    # assert value
    r0 = np.array([0,0,0])
    r_test = pos['tr'][8].r_to_loc(r0)
    rdiff = asnumpy(r_test) - np.array([1,0,-3])
    qry = np.sum(rdiff**2)
    assert(qry<_eps)
    print('passed') 
# test linear scan
    print('Linear scan ... '.format(),end='')    
    #step_size = 0.5
    origin = [0,-1,-1]
    axis = [1,-0.5, -1]
    a0 = np.linalg.norm(axis)
    step_size = 0.5*a0
    steps = list(np.linspace(-3,3,num=7))
    scan = ScanGeometry.linear(steps, origin=origin, axis=axis, 
                               step_size=step_size)
    scan.transform = {'orig':[0,0,0], 'angles':[0, 90, 0], 'order':1}
    pos = scan.positions()
    # scan points
    r = asnumpy(pos['r'])
    
    # ki, kf vectors, assume ki || z, kf || x
    x_ax = np.array([1,0,0])
    z_ax = np.array([0,0,1])
    a = {'ki':[], 'kf':[]}
    for i in range(scan.nsteps):
        x_axt = asnumpy(pos['tr'][i].v_to_loc(x_ax))
        z_axt = asnumpy(pos['tr'][i].v_to_loc(z_ax))
        kf = np.array([r[:,i], r[:,i]+x_axt]).T
        ki = np.array([r[:,i]-z_axt, r[:,i]]).T
        a['kf'].append(kf)
        a['ki'].append(ki)    
    
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
       #ax[i].plot(r[ix,0], r[iy,0], color='red', marker='o', 
       #            linestyle='none', markersize=5) 
        
        ax[i].plot(origin[ix], origin[iy], color='red', marker='o', linestyle='none', 
                markersize=5)
        
        for j in range(scan.nsteps):
            ax[i].plot(a['kf'][j][ix], a['kf'][j][iy], color='red', linestyle='-', 
                    linewidth=0.5)
            ax[i].plot(a['ki'][j][ix], a['ki'][j][iy], color='blue', linestyle='-', 
                    linewidth=0.5)
        ax[i].set_xlabel(lblx)
        ax[i].set_ylabel(lbly)
        ax[i].grid()
        #ax[i].invert_xaxis()
    title = 'origin: '+3*'{:g} '+'axis: '+3*'{:g} '+'step: {:g}'
    fig.suptitle(title.format(*origin, *axis, step_size/a0), fontsize=8)
    fig.tight_layout(w_pad=0.5)
    plt.show() 
    # assert value
    r0 = np.array([0,0,0])
    r_test = pos['tr'][4].r_to_loc(r0)
    rdiff = asnumpy(r_test) - np.array([0.5,-1.25, -1.5])
    qry = np.sum(rdiff**2)
    assert(qry<_eps)
    print('passed') 


class TestObj():
    """Generic class on which the tests are carried out.
    
    Define the surface geometry, test events and test aswers.
    
    This base class creates a "moon" shape from 2 cylinders projected along y.
    
    """
    def __init__(self, title='Moon', scale=1, limits=None, **kwargs):
        # default options
        self.options = {'nev':20, 'proj':[0,2], 'random': True}
        # update from arguments
        for o in self.options:
            if o in kwargs:
                self.options[o] = kwargs[o]
        self.title = title
        self.scale = scale
        self.root = self.create_group()
        self.events = [None, None, None] # [r, vi, vf] events
        self.rays = [] # list of rays for testing
        self._limits = limits

    def _add_surface_limits(self, el, xlim, ylim):
        p = self.options['proj']
        if isinstance(el, Cylinder):
            o = asnumpy(el.trg.orig[p,0])
            xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
            xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
            ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
            ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        elif isinstance(el, ECylinder):
            o = asnumpy(el.trg.orig[p,0])
            xlim[0] = min(xlim[0],1.05*(-el.a+o[0]))
            xlim[1] = max(xlim[1],1.05*(el.a+o[0]))
            ylim[0] = min(ylim[0],1.05*(-el.b+o[1]))
            ylim[1] = max(ylim[1],1.05*(el.b+o[1]))
        elif isinstance(el, Sphere):
            o = asnumpy(el.trg.orig[p,0])
            xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
            xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
            ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
            ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        elif isinstance(el, Plane):
            n = asnumpy(el.trg.v_to_glob(el.n))[p]
            xlim[0] = min(xlim[0],1.05*el.d*n[0])
            xlim[1] = max(xlim[1],1.05*el.d*n[0])
            ylim[0] = min(ylim[0],1.05*el.d*n[1])
            ylim[1] = max(ylim[1],1.05*el.d*n[1])
        return xlim, ylim
    
    def _add_group_limits(self, group, xlim, ylim):
        for el in group.surfaces:
            if isinstance(el, Group):
                xlim, ylim = self._add_group_limits(el, xlim, ylim)
            else:
                xlim, ylim = self._add_surface_limits(el, xlim, ylim)
        return xlim, ylim         
    

    def _create_random_events(self):
        n = self.options['nev']
        p = self.options['proj']
        xlim, ylim = self.limits
        x0 = 0.5*(xlim[1]+xlim[0])
        dx = xlim[1]-xlim[0]
        y0 = 0.5*(ylim[1]+ylim[0])
        dy = ylim[1]-ylim[0]
    
        rn = cp.random.random((3,n)) - 0.5
        r = [0,0,0]
        vi = [0,0,0]
        vf = [0,0,0]
        for i in range(3):
            if i==p[0]:
                r[i] = dx*rn[p[0],:] + x0
                vf[i] = cp.ones(n)
                vi[i] = cp.zeros(n)
            elif i==p[1]:
                r[i] = dy*rn[p[1],:] + y0
                vi[i] = cp.ones(n)
                vf[i] = cp.zeros(n)
            else:
                r[i] = cp.zeros(n)
                vi[i] = cp.zeros(n)
                vf[i] = cp.zeros(n)
        r = cp.array(r)
        # +- 10% random spread of velocities
        dv = 1+(cp.random.random() - 0.5)*0.2
        vi = dv*cp.array(vi)
        vf = dv*cp.array(vf)
        return r,vi,vf    

    def _create_test_events(self):
        
        r =  [[0.5,0,0], [-0.25,0,0.75]]
        vi = [[1,0,2], [1,0,1]]
        vf = [[1,0,0], [1,0,-1]]

        # normalize velocities
        for i in range(len(vi)):
            vn = np.linalg.norm(vi[i])/self.scale
            vi[i] = cp.array(vi[i])/vn
            vn = np.linalg.norm(vf[i])/self.scale
            vf[i] = cp.array(vf[i])/vn
            r[i] = cp.array(r[i])*self.scale
        
        r = cp.asarray(r).T
        vi = cp.asarray(vi).T
        vf = cp.asarray(vf).T
        return r,vi,vf

    def _create_events(self):
        if self.options['random']:
            r,vi,vf = self._create_random_events()
        else:
            r,vi,vf = self._create_test_events()
        self.events = [r, vi, vf]

    def _get_rays(self, traces=None):
        """Collect incident and output rays.
        
        Previous call of :meth:`evaluate` is assumed.
        
        Parameters
        ----------
        traces : list
            Optional list of indices for rays to be selected.
        
        Return
        ------
        list
            A list of ray nodes with shape [2,3]. The first index denotes
            inward/outward cross-section, the 2nd index is for the 
            coresponding coordinates.
        
        """
        root = self.root       
        # select subset of events for plotting 
        if traces is None:
            r, vi, vf = self.events
            is_in = root.is_in
        else:
            r = self.events[0][:,traces]
            vi = self.events[1][:,traces]
            vf = self.events[2][:,traces]
            is_in = root.is_in[traces]       
        
        # collect all incident rays
        out = {'in':[], 'out':[]}
        for path in ['in','out']:
            tx, mx = root._cross_times()[path].values()
            tx = cp.array(tx)
            mx = cp.array(mx)
            rays = []
            if traces is None:
                _min,_mout = mx
                _tin,_tout = tx
            else:
                _min = mx[0][:,traces]
                _tin = tx[0][:,traces]
                _mout = mx[1][:,traces]
                _tout = tx[1][:,traces]               
            if path=='in':
                q_in = cp.array(_tin<0, dtype=int)*_min
                q_out = cp.array(_tout<0, dtype=int)*_mout
            else:
                q_in = cp.array(_tin>0, dtype=int)*_min
                q_out = cp.array(_tout>0, dtype=int)*_mout
            for i in range(r.shape[1]):
                for j in range(len(_tin)):
                    if is_in[i]:
                        if path=='in':
                            if q_in[j,i]:
                                rin = (r[:,i] + vi[:,i]*_tin[j,i])
                                if q_out[j,i]:
                                    rout = (r[:,i] + vi[:,i]*_tout[j,i])
                                else:
                                    rout = r[:,i]
                                rays.append([rin,rout])
                        else:
                            if q_out[j,i]:
                                if q_in[j,i]:
                                    rin = (r[:,i] + vf[:,i]*_tin[j,i])
                                else:
                                    rin = r[:,i]
                                rout = (r[:,i] + vf[:,i]*_tout[j,i])
                                rays.append([rin,rout])
            out[path] = rays
        return out


    def _get_limits(self):
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

    @property
    def limits(self):
        """Limits of view area for plotting and events generation."""
        if self._limits is None:
            return self._get_limits()
        else:
            return self._limits

    def assign_events(self, events=None):
        if events is None:
            self._create_events()
        else:
            self.events = events

    def initiate(self, traces=None):
        """Evaluate events for testing."""
        #self.root.set_events(*self.events) 
        self.root.evaluate(*self.events)
        self.rays = self._get_rays(traces=traces)

    def create_group(self):
        """Create root surface group."""
        sf1 = {'shape':'Cylinder', 'R':self.scale*1, 'axis':'y', 'inner':1,
               'tr':{'orig':[-1.0*self.scale, 0.0, 0.0]}}           
        sf2 = {'shape':'Cylinder', 'R':self.scale*1, 'axis':'y'}
        
        sf1 = {'shape':'ECylinder', 'a':2, 'b':1, 'axis':'y', 'inner':-1,
               'tr':{'orig':[0, 0.0, 0.0]}} 
        root = RootGroup(surfaces=[sf1, sf2], op='and')
        
        return root

    def test_rays(self, n=1000):
        """Test that all points on the ray lines are inside the root group.
        
        The test is done by generating n points on each ray line. For each 
        point, test that it is inside the root surface group.
        
        Return
        ------
        dict
            Two elements 'in', 'out' record results for input/output rays.    
            Each value contains a list of int arrays with the results of the 
            'inside' test (0|1).
        
        """
        eps = 1e-2
        out = {'in':[], 'out':[]}
        t = cp.linspace(eps,1.0-eps, num=10)
        for path in ['in','out']:
            for r in self.rays[path]:
                d = cp.linalg.norm(r[1] - r[0])
                # skip test for trajectories < 0.1 um
                if d>1e-4:
                    v = r[0].reshape(3,1) + t*(r[1] - r[0]).reshape(3,1)
                    #a = cp.asarray(v).T
                    m = self.root.is_inside(r=v)
                    out[path].append(m)
        return out


class TestSimple(TestObj):
    def __init__(self, param, **kwargs):
        self.param = param
        super().__init__(**kwargs)
    
    def create_group(self):
        return RootGroup(**self.param)
        #return Box(size=[1,1,1])
    

class TestCone(TestObj):
    """Test cone shape.
    
    The testing shape is a double-cone connected with a cylinder.
        
    """
    
    def __init__(self, title='Cone', R=0.5, a=0.0, angle=30, **kwargs):
        self.R = R
        self.a = a
        self.angle = angle
        super().__init__(title=title, **kwargs)                
        

    def _get_limits(self):
        shx = self.a*self.scale*np.cos(self.angle*np.pi/180)
        shy = -self.a*self.scale*np.sin(self.angle*np.pi/180)
        ymax = self.R*self.scale
        rho = np.sqrt(ymax**2 + self.scale**2)*1.05
        xlim = [shx-rho, shx+rho]
        ylim = [shy-rho, shy+rho]
        return xlim, ylim

    def _create_test_events(self):
        u = cp.array([[1,0,0],[0,1,0],[0,0,1]])
        p = self.options['proj']
        sg = [-1,0,1]
        events = []
        for i in sg:
            for j in sg:
                r = u[p[0]]*i + u[p[1]]*j
                #r = cp.array([i,j,k])
                events.append(0.25*self.scale*r)
        
        #DEBUG
        #events = [cp.array([0,0,0])]
        ne = len(events)
        
        
        vi = cp.asarray(ne*[u[p[1]]]).T
        vf = cp.asarray(ne*[u[p[0]]]).T
        r = cp.asarray(events).T
        return r,vi,vf

    def create_group(self, inner=-1):
        """Double-cone connected with a cylinder."""
        sc = self.scale
        px1 = {'shape':'Plane', 'inner':-1, 'p':[1,0,0],'d':(1+self.a)*sc}  
        px2 = {'shape':'Plane', 'inner':1, 'p':[1,0,0],'d':(-1+self.a)*sc}    
        px3 = {'shape':'Plane', 'inner':-1, 'p':[1,0,0],'d':(0.5+self.a)*sc}  
        px4 = {'shape':'Plane', 'inner':1, 'p':[1,0,0],'d':(-0.5+self.a)*sc}
        cyl = {'shape':'Cylinder', 'R':0.25*self.R*sc, 'axis':'x', 'inner':-1}        
        tr = {'orig':[0, 0, 0], 'angles':[0,0,0]}
        cone = {'shape':'Cone', 'a':self.a*sc, 'R':self.R, 'axis':'x', 
                'inner':-1, 'tr':tr}
        g_cyl = {'surfaces':[px3,px4, cyl], 'op':'and', 'inner':-1}
        g_cone = {'surfaces':[px1,px2, cone], 'op':'and', 'inner':-1}
        root = RootGroup(surfaces=[g_cyl,g_cone], op='or', 
                             tr={'angles':[0,self.angle,0]})
        return root

    
class TestHexagon(TestObj):
    """Test hexagon composed of cylinders or spheres in a box frame."""
    
    def __init__(self, title='Hexagon', param={}, shape='Cylinder',
                 frame='circle', **kwargs):
        
        self.param = {'R1':1.5, 'R2':2, 'op':['or','or','or','and', 'and'],
                      'inner': [1,1,1,-1, -1], 'rho':2}
        
        for p in self.param:
            if p in param:
                self.param[p] = param[p]
        self.frame = frame
        self.shape = shape
        title = '{}: {}'.format(title, self._get_title())
        super().__init__(title=title, **kwargs)

    def _get_title(self):
        op = self.param['op']
        inner = self.param['inner']
        sgn = ['-','','']
        s = '{}[ '.format(sgn[inner[3]+1])
        for i in range(3):
            s += '{}(S{} {} S{})'.format(sgn[inner[i]+1], 2*i+1, op[i], 2*i+2)
            if i<2:
                s += ' {} '.format(op[3])
        s += ' ] {} {}S7'.format(op[4], sgn[inner[4]+1])
        return s        

    def _create_test_events(self):
        u = cp.array([[1,0,0],[0,1,0],[0,0,1]])
        sg = [-1,0,1]
        lst = []
        for i in sg:
            for j in sg:
                for k in sg:
                    r = cp.array([i,j,k])
                    lst.append(r)
        #DEBUG
        #lst = [cp.array([-1,0,0])]
        group = self.root.surfaces[0]
        events = []
        for g in group.surfaces:
            for s in g.surfaces:
                for p in lst:
                    ev = s.trg.r_to_glob(0.75*p)
                    events.append(ev)
        
        ne = len(events)
        p = self.options['proj']
        vi = cp.asarray(ne*[u[p[1]]]).T
        vf = cp.asarray(ne*[u[p[0]]]).T
        r = cp.asarray(events).T
        return r,vi,vf
    
    def create_group(self):
        """Combine 6 surfaces in 3 groups + outer surface."""
        R1, R2, op, inner, rho = self.param.values()
        # hexagon parameter
        co = np.cos(30*np.pi/180)
        si = 0.5
        # a pair of spheres
        if self.shape=='Sphere':
            sf1 = [{'shape':'Sphere','R':R1,'tr':{'orig':[0.0, 0.0, 0.0]}},
                   {'shape':'Sphere','R':R1,'tr':{'orig':[0.0, 0.0, 2*rho*si]}}
                  ]
        elif self.shape=='Box':
            size = [2*R1, 2*R1, 2*R1]
            ang1 = [0,0,0]
            ang2 = [0,60,0]
            sf1 = [{'shape':'Box','size':size,  
                    'tr':{'orig':[0.0, 0.0, 0.0], 'angles':ang1}},
                   {'shape':'Box','size':size,
                    'tr':{'orig':[0.0, 0.0, 2*rho*si], 'angles':ang2}}
                  ]
        elif self.shape=='ECylinder':
            sf1 = [{'shape':'ECylinder', 'axis':'y', 'a':R1, 'b':R1*0.75, 
                    'angle':-30, 'tr':{'orig':[0.0, 0.0, 0.0]}},
                   {'shape':'ECylinder', 'axis':'y', 'a':R1, 'b':R1*0.75, 
                    'angle':30,'tr':{'orig':[0.0, 0.0, 2*rho*si]}}
                  ]
        
        elif self.shape=='PolygonBar':
            ang1 = [0,0,0]
            ang2 = [0,0,0]
            side=1.0
            height=50.0
            sf1 = [{'shape':'PolygonBar',
                    'num':6, 'side':side,'height':height, 'axis':'y',
                    'tr':{'orig':[0.0, 0.0, 0.0], 'angles':ang1}},
                   {'shape':'PolygonBar',
                    'num':5, 'side':side,'height':height, 'axis':'y',
                    'tr':{'orig':[0.0, 0.0, 2*rho*si], 'angles':ang1}}
                  ]
        else:
            sf1 = [{'shape':'Cylinder', 'axis':'y', 'R':R1,
                    'tr':{'orig':[0.0, 0.0, 0.0]}},
                   {'shape':'Cylinder', 'axis':'y', 'R':R1, 
                    'tr':{'orig':[0.0, 0.0, 2*rho*si]}}
                  ]
        gr = 3*[0]
        # hexagon from 3 pairs of spheres
        gr[0] = {'surfaces':sf1, 'op':op[0], 'inner':inner[0], 
                 'tr':{'orig':[rho*co, 0., -rho*si]}, 
                 'color':(0.5, 0.0, 0.0, 0.15) }   
        
        gr[1] = {'surfaces':sf1, 'op':op[1], 'inner':inner[1], 
                 'tr':{'orig':[0, 0., rho], 'angles':[0, -120, 0]},
                 'color':(0.0, 0.5, 0.0, 0.15)}

        gr[2] = {'surfaces':sf1, 'op':op[2], 'inner':inner[2], 
                 'tr':{'orig':[-rho*co, 0., -rho*si], 'angles':[0, -240, 0]},
                 'color':(0.0, 0.0, 0.5, 0.15)}
        
        
        # we could just rotate around zero origin, but this will test the full
        # transformation
        gr2 = {'surfaces':gr,'op':op[3], 'inner':inner[3], 
                'tr':{'order':1, 'orig':[2.0, 0.0, -2.0], 'rotctr':[-2,0,0], 
                      'angles':[0, -90, 0]}
                }

        if self.frame=='circle':
            # sphere through hexagon corners
            if self.shape=='Sphere':
                sf2 = {'shape':'Sphere', 'R':R2, 'inner':inner[4]}
            else:
                sf2 = {'shape':'Cylinder', 'axis':'y', 'R':R2, 
                        'inner':inner[4]}
            root = RootGroup(surfaces=[gr2,sf2], op=op[4])
        else:
            g1 = {'shape':'Box', 'size':[2*R2,2*R2,2*R2], 'inner':inner[4]}
            g2 = {'shape':'Box', 'size':[2*R2,2*R2,2*R2], 'inner':inner[4],
                  'tr':{'angles':[0,45,0]}}
            # octagone box 
            #px1 = {'shape':'Plane', 'inner':inner[4], 'p':[1,0,0],'d':R2}  
            #px2 = {'shape':'Plane', 'inner':-inner[4], 'p':[1,0,0],'d':-R2}  
            #pz1 = {'shape':'Plane', 'inner':inner[4], 'p':[0,0,1],'d':R2}  
            #pz2 = {'shape':'Plane', 'inner':-inner[4], 'p':[0,0,1],'d':-R2}           
            #g1 = {'surfaces':[px1,px2,pz1,pz2], 'op':'and', 'inner':inner[4]}
            #g2 = {'surfaces':[px1,px2,pz1,pz2], 'op':'and', 'inner':inner[4],
            #      'tr':{'angles':[0,45,0]}}
            root = RootGroup(surfaces=[gr2, g1, g2], op=op[4], 
                         tr={'angles':[0,30,0]})
        return root


  
class TestWires(TestObj):
    """Test holes/wires in a tube."""
    
    def __init__(self, title='Wires', R1=16, R2=3, rho=11, comp='rad', **kwargs):
        
        self.R1 = R1
        self.R2 = R2
        self.comp = comp
        self.rho = rho
        super().__init__(title=title, **kwargs)     

    def _create_test_events(self):
        u = cp.array([[1,0,0],[0,1,0],[0,0,1]])
        sg = [-1,0,1]
        lst = []
        for i in sg:
            for j in sg:
                for k in sg:
                    r = cp.array([i,j,k])
                    lst.append(r)
        group = self.root.surfaces[0]
        events = []
        for g in group.surfaces:
            for s in g.surfaces:
                for p in lst:
                    ev = s.trg.r_to_glob(0.7*p)
                    events.append(ev)
        ne = len(events)
        p = self.options['proj']
        vi = cp.asarray(ne*[u[p[1]]]).T
        vf = cp.asarray(ne*[u[p[0]]]).T
        r = cp.asarray(events).T
        return r,vi,vf
    
    def create_group(self):
        """Three cylindric holes in a rod sample."""
        # define cell
        # angle between rods
        dom = 60
        co = np.cos(dom*np.pi/180)
        si = np.sin(dom*np.pi/180)
        # radial position of rods
        self.rho = 11*self.R1/16

        sf1 = {'shape':'Cylinder', 'axis':'y', 'R':self.R1, 'inner':-1, 
               'tr':{'orig':[0.0, 0.0, 0.0]}}
        sf2 = {'shape':'Cylinder', 'axis':'y', 'R':self.R2, 'inner':1, 
               'tr':{'orig':[self.rho, 0.0, 0.0]}}
        sf3 = {'shape':'Cylinder', 'axis':'y', 'R':self.R2, 'inner':1, 
               'tr':{'orig':[self.rho*co, 0.0, self.rho*si],'angles':[0,dom,0]}}
        sf4 = {'shape':'Cylinder', 'axis':'y', 'R':self.R2, 'inner':1, 
               'tr':{'orig':[self.rho*co, 0.0, -self.rho*si],'angles':[0,-dom,0]}}
        root = RootGroup(surfaces=[sf1, sf2, sf3, sf4], op='and')
        return root
    
    
class TestIntegral():
    """Test convolution integrals on given surface shape.
    
       The surface and other properties are defined by an instance of TestObj.
       
       Parameters
       ----------
       obj : TestObj
           Test object which defines the surface shape and other properties. 
       scan : ScanGeometry
           Scan definition object.
       steps : list
           List of step indices at which teh events are shown. If not given, 
           assume 5 positions evenly distributed along the scan.
           
       Options
       -------
       npos : int
            If steps not given in the constructor, show events at `npos`
            positions evenly distributed along the scan.
       nev_show : int
            Number of events to show at each step. 
       random : bool
           Show events randomly selected from the sampling file, otherwise
           select required number of events sequentially from the beginning.
    
    """
    
    def __init__(self, obj, scan, steps=None, reference=None, efunc=None, 
                 **kwargs):
        self.obj = obj
        self.scan = scan
        self.root = obj.root
        self.steps = steps
        self.reference = reference
        self.options = {'random': False, 'nev':5000, 'npos':5, 'nev_show':20}
        for key in kwargs:
            if key in self.options:
                self.options[key] = kwargs[key]
            
        self._efunc = efunc

    def _get_traces(self):
        """Generate scan positions in step units for plotting of traces.
        
        Parameters
        ----------

        """
        n = self.options['nev_show']
        idx = np.linspace(0, n, num=n, endpoint=False, dtype=int)
        if self.steps is None:
            nsteps = self.scan.nsteps
            npos = self.options['npos']
            self.steps = [(1+int(nsteps/npos))*i for i in range(npos)]
        ns = len(self.steps)
        res = np.zeros(n*ns, dtype=int)
        for i in range(ns):
            res[i*n:(i+1)*n] = idx + self.steps[i]*self.integ.nev
        return res                


    def initiate(self, sampling, material):
        self.sampling = sampling
        self.material = material
        self.cell = Cell(self.root, self.material)
        self.cell.reference = self.reference
        # define integral and events x scan steps 
        self.integ = GaugeIntegral(self.cell, self.sampling, self.scan)
        # NOTE - this calls set_events and evaluate on obj.root 
        self.integ.initiate(nev=self.options['nev'], 
                            random=self.options['random'])
        # set integration events (replaces _create_events)
        self.obj.assign_events(events=[self.integ._ev['r'], 
                                    self.integ._ev['ki'], 
                                    self.integ._ev['kf']])       

    def efunc(self, x):
        """Intrinsic strain function."""
        if self._efunc is None:
            return 0
        else:
            return self._efunc(x)

    def run_test(self):
        # make convolution integrals 
        gauge = self.integ.conv_gauge()
        conv = self.integ.conv_strain(efunc=self.efunc)
        return to_numpy(gauge), to_numpy(conv)

    
#%% Test performers
class AbstractTest():    
    
    def __init__(self):
        self.root = None # root surface group  
        self.traces = None # optional list of indices for trace plot   
        self.tests = {} # list of test objects
        self.obj = None # actually runing TestObj object
        self.options = {}
        
    def _add_surface_limits(self, el, xlim, ylim):
        if isinstance(el, Cylinder):
            o = asnumpy(el.trg.orig[0::2,0])
            xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
            xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
            ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
            ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        if isinstance(el, Sphere):
            o = asnumpy(el.trg.orig[0::2,0])
            xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
            xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
            ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
            ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        elif isinstance(el, Plane):
            n = asnumpy(el.trg.v_to_glob(el.n))
            xlim[0] = min(xlim[0],1.05*el.d*n[0])
            xlim[1] = max(xlim[1],1.05*el.d*n[0])
            ylim[0] = min(ylim[0],1.05*el.d*n[2])
            ylim[1] = max(ylim[1],1.05*el.d*n[2])
        return xlim, ylim

    def _update_options(self, **kwargs):
        for o in self.options:
            if o in kwargs:
                self.options[o] = kwargs[o]

    
    def _add_group_limits(self, group, xlim, ylim):
        # assume (x,z) projection plane
        for el in group.surfaces:
            if isinstance(el, Group):
                xlim, ylim = self._add_group_limits(el, xlim, ylim)
            else:
                xlim, ylim = self._add_surface_limits(el, xlim, ylim)
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

    def plot_cell(self, ax, limits=[[-1,1],[-1,1]], npix=(500,500)):
        xlim, ylim = limits
        mask = self.obj.root.get_map(xlim=xlim, ylim=ylim, npix=npix)
        cmap = mpl.colormaps['binary']
        extent = np.append(xlim, ylim)
        ax.imshow(mask, cmap=cmap, alpha=0.5, vmin=0, vmax=3,
                  origin='lower', extent=extent)

    def plot_trace(self, ax):
        """Plot incident and output rays."""
        lines = {'in':'b-','out':'r-'}
        # collect all incident rays
        for path in ['in','out']:
            if len(self.obj.rays[path])>0:
                rays = asnumpy(cp.asarray(self.obj.rays[path]))
                cross_in = rays[:,0,:]
                cross_out = rays[:,1,:]
                ax.plot(cross_in[:,0], cross_in[:,2], 'b.', markersize=2)
                ax.plot(cross_out[:,0], cross_out[:,2], 'r.', markersize=2)
                for i in range(rays.shape[0]):
                    x = rays[i,:,0]
                    y = rays[i,:,2]
                    ax.plot(x, y, lines[path], linewidth=0.5)


    def plot(self, title='', trace=True, limits=None):
        fnt = 12
        if limits:
            xlim, ylim = limits
        else:
            xlim, ylim = self.get_limits()
        fig, ax = plt.subplots(figsize=(7,7))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.plot_cell(ax, limits=[xlim, ylim])
        if self.traces is None:
            r = self.obj.events[0]
            is_in = self.obj.root.is_in
        else:
            r = self.obj.events[0][:,self.traces]
            is_in = self.obj.root.is_in[self.traces]
        r0 = asnumpy(r)
        if trace:
            m = asnumpy(is_in)
        else:
            m = np.ones(r0.shape[1])
        r1 = r0[:,(m==0)]
        r2 = r0[:,(m>0)]
        ax.plot(r1[0], r1[2], color='gray', marker='x', linestyle='none', markersize=3)
        ax.plot(r2[0], r2[2], color='black', marker='x', linestyle='none', markersize=4)
        if trace:
            self.plot_trace(ax)
        ax.set_xlabel('x',fontsize=fnt)
        ax.set_ylabel('z',fontsize=fnt)
        ax.grid()
        if title:
            plt.title(title, fontsize=fnt)
        plt.show() 
    
    def run(self, tests='all', **options):
        """Run requested tests.
        
        Parameters
        ----------
        tests : str or list
            Which test to run, either 'all' or a list. The list 
            can contain aither test names or indexes from the pre-defined list.
        **options : 
            Test options.
        """       
        self._update_options(**options)
        self.define_tests(**options)
        if isinstance(tests,list):
            torun = {}
            keys = list(self.tests.keys())
            for i in tests:
                if isinstance(i,int):
                    key = keys[i]
                else:
                    key = i
                torun[key] = self.tests[key]
        else:
            torun = self.tests        
        for key in torun:
            print('Run {} '.format(key), end='')
            try:
                self.run_test(key)
                print('passed')
                # check also nvidia-smi
                mempool = cp.get_default_memory_pool()
                mb = 1024**2
                print('memory used: {:g}'.format(mempool.used_bytes()/mb))
                print('memory total: {:g}'.format(mempool.total_bytes()/mb))
            except Exception as e:
                print('failed')
                mempool = cp.get_default_memory_pool()
                mb = 1024**2
                print('memory used: {:g}'.format(mempool.used_bytes()/mb))
                print('memory total: {:g}'.format(mempool.total_bytes()/mb))
                raise(e)   
        
    @abc.abstractmethod
    def run_test(self):
        """Define test operations."""
        msg = 'Methot run_test() not defined on AbstractText class.'
        raise Exception(msg)
    
class TestSurface(AbstractTest):
    def __init__(self):
        super().__init__()
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
        self.options.update({'composite':False, 'plot':True, 'nev':20, 
                        'random':False})
    
    def _plot_surface_colored(self, ax, surf, color):  
        # assume Cylinder || y in (x,z) plane
        o = asnumpy(surf.trg.orig[0::2,0])
        x = surf.R*self.si + o[0]
        y = surf.R*self.co
        y2 = y + o[1]
        y1 = -y + o[1]
        ax.fill_between(x, y1, y2, facecolor=color)   
        
    def _plot_group_colored(self, ax, group, color):
        for el in group.surfaces:
            if isinstance(el, Group):
                cl = color
                if hasattr(el,'color'):
                    cl = el.color
                self._plot_group_colored(ax, el, cl)
            else: 
                self._plot_surface_colored(ax, el, color)
        


    def plot_cell(self, ax, **kwargs):
        """Plot cell as a composite of surface groups, or a pixel map."""
        if self.options['composite']:
            gray = (0.5, 0.5, 0.5, 0.15)
            self._plot_group_colored(ax, self.root, gray)
        else:
            super().plot_cell(ax, **kwargs)

    def define_tests(self, **options):
        self.tests.clear()  
        
        orig = [3,0,0]
        rotctr = [-3,0,0]
        stp = -60
        frame = {'shape':'Sphere', 'R':4, 'inner':-1} 
        
        cyl = {'shape':'Cylinder', 'R':1, 'axis':'y', 'inner':1, 
               'tr':{'orig':orig, 'rotctr':rotctr, 'angles':[0,4*stp,0]}} 
        sph = {'shape':'Sphere', 'R':1, 'inner':1, 
               'tr':{'orig':orig, 'rotctr':rotctr, 'angles':[0,1*stp,0]}} 
        box = {'shape':'Box', 'size':[1.75,1.25,1.5], 'inner':1,
               'tr':{'orig':orig, 'rotctr':rotctr, 'angles':[0,2*stp,0]}} 
        
        cone = {'shape':'Cone', 'a':1, 'axis':'y', 'inner':1, 
               'tr':{'orig':[3,-1, 0], 'rotctr':rotctr, 'angles':[0,3*stp, 30]}} 

        erod = {'shape':'ECylinder', 'a':1, 'b':0.5, 'axis':'y', 'inner':1, 
                'angle':90,
               'tr':{'orig':orig, 'rotctr':rotctr, 'angles':[0,5*stp,0]}}
        pbar = {'shape':'PolygonBar', 'num':6, 'side':1, 'axis':'y', 'inner':1,
               'tr':{'orig':orig, 'rotctr':rotctr, 'angles':[0,0,0]}}

        g = {'surfaces':[cyl,sph,box,cone, erod, pbar],'op':'and'}
        param = {'surfaces':[frame, g], 'op':'and'}

        self.tests['Primitives'] = TestSimple(title='Primitives', param=param, 
                                              **options)
        
        
        self.tests['Moon'] = TestObj(title='Moon', **options)
        
        param = {'R1':0.8, 'R2':2, 
                 'op':['or','or','or','and','and'], 
                 'inner':[1,1,1,-1,-1]}
        self.tests['Hexagon1'] = TestHexagon(title='Hexagon1', param=param, 
                                             shape='Sphere', **options)
        
        #DEBUG 'R1':1.2
        param = {'R1':1.2, 'R2':2, 
                 'op':['or','or','or', 'or','and'], 
                 'inner':[-1,-1,-1,-1, 1]}  
        self.tests['Hexagon2'] = TestHexagon(title='Hexagon2', param=param, 
                                             **options)
        param = {'R1':1.2, 'R2':2, 
                 'op':['or','or','or', 'or','and'], 
                 'inner':[-1,-1,-1,-1, -1]}  
        self.tests['Hexagon3'] = TestHexagon(title='Hexagon3', param=param, 
                                             frame='box', **options)

        param = {'R1':1.2, 'R2':2, 
                 'op':['and','or','and', 'or','and'], 
                 'inner':[-1,-1,-1,-1, -1]}  
        self.tests['Hexagon4'] = TestHexagon(title='Hexagon4', param=param, 
                                             **options)

        param = {'R1':0.5, 'R2':3, 
                 'op':['or','or','or', 'or','and'], 
                 'inner':[-1,-1,-1, 1,-1]}  
        self.tests['Hexagon_Box'] = TestHexagon(title='Hexagon_Box', 
                                                param=param, shape='Box', 
                                                limits = [[-4,4],[-4,4]],
                                                **options)
        param = {'R1':0.5, 'R2':3, 
                 'op':['or','or','or', 'or','and'], 
                 'inner':[-1,-1,-1, 1,-1]}  
        self.tests['Hexagon_Ell'] = TestHexagon(title='Hexagon_Ell', 
                                                param=param, shape='ECylinder', 
                                                limits = [[-4,4],[-4,4]],
                                                **options)

        self.tests['Cone'] = TestCone(title='Cone', angle=30, a=-0.5, scale=2, 
                                      **options)

    def run_test(self, name): 
        self.obj = self.tests[name]
        self.obj.assign_events()
        self.obj.initiate()      
        #self.root = obj.root     
        #self.r = obj.events[0]
        #self.vi = obj.events[1]
        #self.vf = obj.events[2]
        if self.options['plot']:
            self.plot(title='test surface: {}'.format(self.obj.title), 
                      limits=self.obj.limits)
        
        q = self.obj.test_rays()
        # for in,out paths
        for p in q:
            # loop for rays
            for i in range(len(q[p])):
                try:
                    assert(all(q[p][i]))
                except Exception as e:
                    rin = self.obj.rays[p][i][0] 
                    rout = self.obj.rays[p][i][1] 
                    delta = rout-rin
                    print('path:    {}'.format(p))
                    print('traj No: {}'.format(i))
                    print('delta: {}'.format(delta))
                    print('rin: {}'.format(rin))
                    print('rout: {}'.format(rout))
                    raise e

    
    def run(self, tests='all', **options):
        """Test tracing through groups of surfaces."""
        print('\nTestSurface:')
        super().run(tests=tests, **options)
        

class TestConv(AbstractTest):
    """Test convolution methods.
    
    See docstrings of :obj:`stressfit.tracing.integrals`.
    
    """
    def __init__(self, save=False):
        super().__init__()
        self.cell = None
        self.save = save
        # TestIntegral object on which the test is actually running
        self.conv_test = None 
        
    def save_gauge(self, gauge, fname='test.csv'):
        if 'cog' in gauge:
            gauge['cogx'] = gauge['cog'][0]
            gauge['cogy'] = gauge['cog'][1]
            gauge['cogz'] = gauge['cog'][2]
            del gauge['cog']
        if 'qref' in gauge:
            gauge['q_X'] = gauge['qref'][0]
            gauge['q_Y'] = gauge['qref'][1]
            gauge['q_Z'] = gauge['qref'][2]
            del gauge['qref']
        df = pd.DataFrame(gauge)
        df.to_csv(fname, index=False)
        
    #override 
    def get_limits(self):
        r = asnumpy(self.obj.events[0])
        xmin = min(r[0,:])
        xmax = max(r[0,:])
        ymin = min(r[2,:])
        ymax = max(r[2,:])
        dr =  [xmax-xmin, ymax-ymin]
        ctr = [0.5*(xmin+xmax), 0.5*(ymin+ymax)]
        scale = max(dr)*1.05
        xlim = [0,0]
        ylim = [0,0]
        xlim[0] = ctr[0] - 0.5*scale 
        xlim[1] = ctr[0] + 0.5*scale 
        ylim[0] = ctr[1] - 0.5*scale 
        ylim[1] = ctr[1] + 0.5*scale 
        return xlim, ylim        


    def define_circular_scan(self, comp='rad', shape='Cylinder', frame='circle'):
        """Test circular scan for rad or hoop component."""   
        param = {'R1':3, 'R2':16, 
                 'op':['or','or','or', 'or','and'], 
                 'inner':[-1,-1,-1, 1,-1],'rho':11}  
        obj = TestHexagon(title='Hexagon_Box', param=param, shape=shape, 
                          frame=frame)
        
        #obj = TestWires(comp=comp)
        # define scan
        nsteps = 51
        steps = np.linspace(-90, 90, num=nsteps)
        #origin_rot = [-a*com,0,a*som]
        origin_rot = [0,0,0]
        origin = [param['rho'],0,0]
        #if comp=='rad':
        #    origin = [-1.2,0,1.2]
        #else:
        #    origin = [1.0,0,-1.0]
        scan = ScanGeometry.circular(steps, origin=origin, 
                                          axis=[0,1,0], 
                                          origin_rot=origin_rot, 
                                          step_size=1.0)        
        # sample rotation
        ang_comp = {'rad':45, 'hoop':-45}
        om = ang_comp[comp] 
        scan.transform = {'angles':[0, om, 0]}
        #steps = [(1+int(nsteps/npos))*i for i in range(npos)]
        steps = [0, int(0.25*nsteps), int(0.5*nsteps), int(0.75*nsteps), nsteps-1]
        def efunc(x):
            eps = 0.6e-3*cp.cos(2*x*cp.pi/180)
            return eps
        
        test = TestIntegral(obj, scan, reference=1, steps=steps, efunc=efunc,
                            nev=10000, nev_show=20)
        return test
        
    def define_tests(self, **options):
        self.tests.clear()
        
        # test 1
        def efunc(x):
            eps = 1e-3*(x-5)/5
            return eps
        obj = TestObj(title='Moon', scale=10, efunc=efunc, **options)
        steps = np.linspace(-1, 11, num=21)
        scan = ScanGeometry.linear(steps, axis=[1,0,0])
        scan.transform = {'angles':[0, -45, 0]}
        test = TestIntegral(obj, scan, reference=0, efunc=efunc)
        self.tests['Test_1'] = test
        
        # test 2
        self.tests['Test_2'] = self.define_circular_scan(comp='rad', 
                                                         shape='Cylinder', 
                                                         frame='circle')

        # test 3
        self.tests['Test_3'] = self.define_circular_scan(comp='hoop',
                                                         shape='Box', 
                                                         frame='circle') 
        
        self.tests['Test_4'] = self.define_circular_scan(comp='hoop',
                                                         shape='ECylinder', 
                                                         frame='circle')    
    def initiate(self):
        # load sampling events
        evdata = load_data('events_S_1mm.dat',kind='data')
        my_data = {'data':evdata,'columns':[1, 4, 7, 10, 11]}
        self.sampling = StrainSampling(my_data)
        self.sampling.print_properties()
        # define material with beam attenuation
        att = load_data('mu_Cu.dat',kind='tables')
        self.material = Material()
        self.material.set_extinction(Extinction.TABLE, table=att)
    
    def errorbar(self, ax, mask, x, y, yerr=None, label=None, **kwargs):
        """Like ax.errorbar, but plot only valid segments."""
        xx = []
        yy = []
        err = []
        i1 = -1
        i2 = -1
        m = mask[0]
        if m:
            i1 = 0
        for i in range(len(x)):
            if m:
                if not mask[i] or i==len(x)-1:
                    i2 = i
                    xx.append(x[i1:i2])
                    yy.append(y[i1:i2])
                    if yerr is not None:
                        err.append(yerr[i1:i2])
                    else:
                        err.append(None)
                    m = False
            if not m:
                if mask[i]:
                    i1 = i
                    m = True
            
        for i in range(len(xx)):
            if i==0:
                ax.errorbar(xx[i], yy[i], yerr=err[i], label=label, **kwargs)
            else:
                ax.errorbar(xx[i], yy[i], yerr=err[i], **kwargs)
        
    
    def plot_conv(self, conv, scan):
        if scan.shape  == ScanGeometry.CIRCULAR:
            ftit = 'Circular scan'
        else:
            ftit = 'Linear scan'
        xtit = 'scan position, {}'.format(scan.units)
        ytit = 'shift of the centre, {}'.format(scan.units)
        
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6))
        xmin = min(conv['x'])
        dx = max(conv['x'])-min(conv['x'])
        xlim = {'left':xmin-0.05*dx, 'right':xmin+1.05*dx}
        
        mask = conv['mask'] # self.integ._gauge['sg']
        # counts
        ax=axs[0,0]
        self.errorbar(ax, mask, conv['x'], conv['cnts'], yerr=conv['err_cnts'], 
                      fmt='ko-', markersize=3)
        ax.set_xlim(**xlim)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit)
        ax.set_ylabel('sampling volume, rel. units')
        ax.grid()
        
        # pseudo-strain
        ax=axs[0,1]
        self.errorbar(ax, mask, conv['x'], 1e6*conv['eps'], yerr=1e6*conv['err_eps'], 
                    fmt='ko-', markersize=3)
        ax.set_xlim(**xlim)
        #ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit)
        ax.set_ylabel('pseudo-strain, 1e-6')
        ax.grid()
        
        # scan position shift and width
        ax=axs[1,0]        
        y = conv['pos'] - conv['x']        
        self.errorbar(ax, mask, conv['x'], y, yerr=conv['err_pos'], fmt='bo-',  
                    markersize=3, label='shift')
        self.errorbar(ax, mask, conv['x'], conv['width'], fmt='rx-',  
                    markersize=3, label='width')
        ax.set_xlim(**xlim)
        ax.legend()
        ax.set_xlabel(xtit)
        ax.set_ylabel(ytit)
        ax.grid()
        
        # cog
        ax=axs[1,1]
        self.errorbar(ax, mask, conv['x'], conv['cog'][0], fmt='ro-', markersize=3, 
                    label='x')
        self.errorbar(ax, mask, conv['x'], conv['cog'][1], fmt='go-', markersize=3, 
                    label='y')
        self.errorbar(ax, mask, conv['x'], conv['cog'][2], fmt='bo-', markersize=3, 
                    label='z')
        if scan.shape == ScanGeometry.CIRCULAR:
            # add radial component for circular scan
            scan_pos = to_numpy(scan.to_scan_coord(conv['cog']))
            self.errorbar(ax, mask, conv['x'], scan_pos['r'], fmt='ko-', 
                          markersize=3, label='r')
        ax.legend()
        ax.set_xlabel(xtit)
        ax.set_ylabel('centre of gravity, mm')
        ax.grid()
        # plot
        fig.suptitle(ftit)
        fig.tight_layout()
        plt.show()

    def plot_strain(self, gauge, strain, scan):
        """Plot strain as a function of (i) scan position and (ii) centre
        of gravity.
        
        Parameters
        ----------
        gauge : dict
            Result of :meth:`GaugeIntegral.conv_gauge`.
        strain : dict
            Result of :meth:`GaugeIntegral.conv_strain`.
        
        """
        if scan.shape  == ScanGeometry.CIRCULAR:
            ftit = 'Circular scan'
        else:
            ftit = 'Linear scan'
        xtit1 = 'scan position, {}'.format(scan.units)
        xtit2 = 'centre of gravity, {}'.format(scan.units)
        ytit = 'strain, 1e-6'
        mask = strain['mask'] #self.integ._gauge['sg']
        fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(4,6))
        
        # x = scan position
        ax = axs[0]
        self.errorbar(ax, mask, strain['x'], 1e+6*strain['eps'], yerr=1e+6*strain['err_eps'], 
                    fmt='ko-', markersize=3, label='total')
        self.errorbar(ax, mask, gauge['x'], 1e+6*gauge['eps'], yerr=1e+6*gauge['err_eps'], 
                    fmt='b--', markersize=3, label='pseudo-strain')
        
        yerr = np.sqrt(gauge['err_eps']**2+strain['err_eps']**2)
        self.errorbar(ax, mask, strain['x'], 1e+6*(strain['eps']-gauge['eps']), yerr=1e+6*yerr, 
                    fmt='r-', markersize=3, label='intrinsic strain')
        xmin = min(gauge['x'])
        dx = max(gauge['x'])-min(gauge['x'])
        xlim = {'left':xmin-0.05*dx, 'right':xmin+1.05*dx}
        ax.set_xlim(**xlim)
        #ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit1)
        ax.set_ylabel(ytit)
        ax.legend()
        ax.grid()        
        # x = scan position
        ax = axs[1]
        self.errorbar(ax, mask, gauge['pos'], 1e+6*strain['eps'], yerr=1e+6*strain['err_eps'], 
                    fmt='ko-', markersize=3, label='total')
        self.errorbar(ax, mask, gauge['pos'], 1e+6*gauge['eps'], yerr=1e+6*gauge['err_eps'], 
                    fmt='b--', markersize=3, label='pseudo-strain')        
        yerr = np.sqrt(gauge['err_eps']**2+strain['err_eps']**2)
        self.errorbar(ax, mask, gauge['pos'], 1e+6*(strain['eps']-gauge['eps']), yerr=1e+6*yerr, 
                    fmt='r-', markersize=3, label='intrinsic strain')
        xmin = min(gauge['pos'])
        dx = max(gauge['pos'])-min(gauge['pos'])
        xlim = {'left':xmin-0.05*dx, 'right':xmin+1.05*dx}        
        ax.set_xlim(**xlim)
        #ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit2)
        ax.set_ylabel(ytit)
        ax.legend()
        ax.grid()        
        # plot
        fig.suptitle(ftit)
        fig.tight_layout()
        plt.show()

    def run_test(self, name):        
        self.conv_test = self.tests[name]
        self.obj = self.conv_test.obj
        self.conv_test.initiate(self.sampling, self.material)
        gauge, conv = self.conv_test.run_test()
        self.traces = self.conv_test._get_traces()
        self.obj.rays = self.obj._get_rays(traces=self.traces)
        self.obj.root.cleanup()
        # plot event map
        self.plot()
        # plot results
        self.plot_conv(gauge, self.conv_test.scan)
        if self.save:
           self.save_gauge(gauge, fname='{}_results.csv'.format(name))
        self.plot_strain(gauge, to_numpy(conv), self.conv_test.scan)
        
    def init_bench(self):
        self.conv_test = self.define_circular_scan(comp='rad', shape='ECylinder')
        self.obj = self.conv_test.obj
        self.conv_test.initiate(self.sampling, self.material)
        
    def bench_func(self):
        #f = np.random.random()
        self.conv_test.integ.conv_gauge()
        self.conv_test.integ.conv_strain(efunc=self.conv_test.efunc)
        #conv = self.conv_test.integ.conv_strain()

    def run(self, tests='all', **options):
        self.initiate()
        print('\nTestConv:')
        super().run(tests=tests, **options)
        
        
    def bench(self, n_repeat=10):
        if has_gpu():
            from cupyx.profiler import benchmark
            #self.initiate()
            print('bench init:')
            self.init_bench()
            print('bench run:')
            res =  benchmark(self.bench_func, (), n_repeat=n_repeat)
            print(res)
            time = cp.average(res.gpu_times) + cp.average(res.cpu_times)
        else:
            import timeit
            glob = globals()
            glob.update({'self':self})
            res = timeit.timeit('self.bench_func()',
                                globals=glob, 
                                setup='print("bench init:");self.init_bench();print("bench run:")',
                                number=n_repeat)
            time = res/n_repeat
        return time

def bench_convolution():
    test = TestConv()
    test.initiate()
    test.init_bench()
    print('GPU: {}'.format(has_gpu()))
    time = test.bench(n_repeat=100)
    print('time: {:g}'.format(time))

#%% Main
if __name__ == "__main__":
    print('test_tracing:')
    for test in _run_tests:
        if test=='primitives':
            test = TestSurface()
            test.run(tests=['Primitives'], nev=100, random=True)
           # test.run(tests=['Hexagon_Ell'], nev=100, random=True)
        if test=='integrals' in _run_tests:
            test = TestConv(save=True)
            test.run(tests=['Test_4'])
        if test=='scan':
            test_scan()
        if test=='bench':
            bench_convolution()
        


#%% Test cone

def test_cone_depth():
    eps = 1e-10
    s2 = 1/np.sqrt(2)
    r = np.array([-s2,-s2,3]) 
    R = 1
    ra = r[2] - 0
    
    # cone angles cos and sin
    co = 1/np.sqrt(1+R**2)
    si = R*co
    # aximuthal component
    x = np.linalg.norm(r[0:2])
    # axial component
    z = np.abs(ra)
    depth = z*si - x*co
    assert(abs(depth-1/s2)<eps)

def test_cone_coord():
    eps = 1e-10
    s2 = 1/np.sqrt(2)
    r = np.array([[s2],[s2],[3]])
    R = 1    
    # cone angles cos and sin
    co = 1/np.sqrt(1+R**2)
    si = R*co
    # aximuthal component
    x = np.linalg.norm(r[0:2])
    
    nr = r.shape[1]
    one = np.ones(nr,dtype=float)
    nul = np.zeros(nr,dtype=float)
    # mask for out-of-axis points
    mx = np.array(x>eps, dtype=int)
    # azimuthal angle of the radial component
    # avoid division-by-zero errors
    # we could use arctan2, but this should be fatser ...
    if np.sum(mx)<nr:
        cosia = mx*r[0:2,:]/(x+(1-mx)*eps) + (1-mx)*cp.array([one, nul])
    else:
        cosia = r[0:2,:]/x
    # sign of z-coordinate
    zsn = np.sign(r[2,:])
    # principal axes for in-plane coordinates (r[1]=0)
      # z = cp.array([co*one,nul,-si*zsn]) # normal to surface
      # y = cp.array([si*one,nul,co*zsn]) # along surface axis
      # x = cp.array([nul,one,nul]) # hoop 
    # rotate by the azimuthal angle around cone axis
    Z = [co*cosia[0,:], co*cosia[1,:], -si*zsn] # normal to surface
    Y = [si*cosia[0,:], si*cosia[1,:], co*zsn]  # along surface axis
    X = [-cosia[1,:], cosia[0,:], nul ]           # hoop 
    U = cp.asarray([X, Y, Z])
    
    assert(abs(U[0,:,0].dot(U[1,:,0])) < eps)
    assert(abs(U[0,:,0].dot(U[2,:,0])) < eps)
    assert(abs(U[1,:,0].dot(U[2,:,0])) < eps)
    x = [-s2, s2, 0]
    y = [0.5, 0.5, s2]
    z = [0.5, 0.5, -s2]
    u = cp.asarray([x,y,z])
    for i in range(3):
        qx = abs(U[i,:,0] - u[i,:])  < eps
        try:
            assert(all(qx))
        except Exception as e:
            print('Error in test {}: {}'.format(i,qx))
            print(U[i,:,0])
            print(u[i,:])
            raise e
    

       
#test_cone_depth()
#test_cone_coord()

#%% Test ellipse

from stressfit.tracing.primitives import ellipse_nearest_point


def draw(a, b, r=None, rc=None, dist=None, plotnorm=True, color='r'):
    def get_normal(rc, a, b):
        mat = np.diagflat([b**2,a**2])
        n = np.dot(mat,rc)
        return n/np.linalg.norm(n, axis=0)
    # plot ellipse
    t = 2*np.pi*np.linspace(0,1,num=100, dtype=float)
    y = b*np.sin(t)
    x = a*np.cos(t)
    fig, ax = plt.subplots(figsize=(7,7))
    rng = max(a,b)
    xlim = [-rng*1.5, rng*1.5]
    ylim = xlim
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.plot(x,y,'k-')
    # plot sample points
    pt_args = {'linestyle':'none', 'marker':'.'}
    rs = None
    if r is not None:
        rs = asnumpy(r)
        d = asnumpy(dist)
        if dist is not None:
            rs_none = rng*10*np.ones(rs.shape)
            m = np.array(d<0, dtype=int)
            rs_in = m*rs + (1-m)*rs_none
            rs_out = m*rs_none + (1-m)*rs
            ax.plot(rs_in[0,:], rs_in[1,:], color='black', **pt_args)
            ax.plot(rs_out[0,:], rs_out[1,:], color='gray', **pt_args)
        else:
            ax.plot(rs[0,:], rs[1,:], color='black', **pt_args)
    
    # plot normals and distances
    line_solid = {'linestyle':'solid', 'linewidth':0.75, 'marker':'none'}
    line_dash = {'linestyle':'solid', 'linewidth':0.5, 'marker':'none'}
    if rc is not None:
        rr = asnumpy(rc)
        ax.plot(rr[0,:], rr[1,:], color=color, **pt_args)
        if rs is not None:
            for i in range(rr.shape[1]):
                ax.plot([rr[0,i],rs[0,i]], [rr[1,i],rs[1,i]], 
                        color=color, **line_solid)
        if plotnorm:
            vn = get_normal(rr, a, b)
            rn = rr + vn 
            for i in range(rr.shape[1]):
                ax.plot([rr[0,i],rn[0,i]], [rr[1,i],rn[1,i]], 
                        color='gray', **line_dash)
    ax.grid()
    plt.show()
 

# TODO: debug, previous working is version in tracing.primitives
def test_ellipse_dist(a=5, b=2, ns=1000, tol=0.001, verbose=False):
    """Test search for nearest point at an ellipse."""
    if verbose:
        msg = 'Test distance from ellipse, a={:g}, b={:g}, tol={:g}, np={:g}'
        print(msg.format(a, b, tol, ns))
    # get points on the ellipse 
    phi = cp.linspace(0, 2*cp.pi, num=ns, endpoint=False)
    tx = cp.cos(phi)
    ty =  cp.sin(phi)
    r0 = cp.array([a*tx,b*ty])
    n = len(phi) # number of points
    nb = max(1,int(ns/100)) # skip points for visualization, show up to 100 pts
    # centers of curvature
    ctr = cp.array([(a*a - b*b)*tx**3/a, 
                    (b*b - a*a)*ty**3/b]) + 1e-10
    # normals to the ellipse
    r1 = r0 - ctr
    # minimum and maximum distances from the ellipse, relative to r1
    kmin = cp.maximum(-ctr[0,:]/r1[0,:], -ctr[1,:]/r1[1,:])
    kmax = cp.ones(n)
    # random points along normals at r0
    cp.random.seed(1001)
    r =  ctr + (kmin + (kmax-kmin)*2*cp.random.rand(n))*r1
    # normals from these points outwards   
    mat = cp.diagflat(cp.array([b**2,a**2]))
    n = cp.dot(mat,r0)
    rn = n/np.linalg.norm(n, axis=0)
    # distance of p along these normals, positive for points outside the ellipse 
    dist0 = cp.sum(rn*(r-r0), axis=0)
    if verbose:
        draw(a, b, r=r[:,::nb], rc=r0[:,::nb], dist=dist0[::nb],color='r')
    # evaluated nearest distance from the ellipse
    dist1, r1, rn1 = ellipse_nearest_point(a, b, r, tol=0.001)
    # dist1, rc1  = ell_nearest_point(a, b, r, verbose=verbose, tol=tol)
    if verbose:
        draw(a, b, r=r[:,::nb], rc=r1[:,::nb], dist=dist1[::nb], color='g')
    # difference in nearest points positions
    dif = r1 - r0
    # difference of distances
    ddist = dist1-dist0
    # check for the maximum difference
    qry1 = cp.max(cp.abs(dif))
    qry2 = cp.max(cp.abs(ddist))
    if verbose:
        print('maximum coordinate difference: {:g}'.format(qry1))
        print('maximum distance difference: {:g}'.format(qry2))
    assert(qry2<tol)
    
#test_ellipse_dist(a=5, b=2, verbose=True, ns=100000, tol=0.001)


