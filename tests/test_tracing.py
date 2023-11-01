# -*- coding: utf-8 -*-
"""
Tests of stressfit.tracing package.

Created on Mon Sep 25 15:33:39 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

#%% Run test
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from stressfit.tracing import options
options['gpu'] = True
from stressfit.tracing.cuda import cp, asnumpy, has_gpu, to_numpy
from stressfit.tracing.cells import Cell, Extinction, Material
from stressfit.tracing.primitives import Transform, Cylinder, Group
from stressfit.tracing.events import StrainSampling
from stressfit.tracing.scan import ScanGeometry
from stressfit.tracing.integrals import GaugeIntegral
from stressfit.dataio import load_data


_run_tests = ['integrals']
#_run_tests = ['bench']

class TracingTest():    
    
    def __init__(self):
        self.root = None # root surface group  
        self.traces = None # optional list of indices for trace plot
        
        
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
        tr = Transform(**kwargs)        
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
                o = asnumpy(el.trg.orig[0::2,0])
                xlim[0] = min(xlim[0],1.05*(-el.R+o[0]))
                xlim[1] = max(xlim[1],1.05*(el.R+o[0]))
                ylim[0] = min(ylim[0],1.05*(-el.R+o[1]))
                ylim[1] = max(ylim[1],1.05*(el.R+o[1]))
        return xlim, ylim        
    
    def set_events(self, r, vi, vf):
            self.r = r
            self.vi = vi
            self.vf = vf
            self.root.set_events(r, vi, vf)  
    
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

    def plot_cell_map(self, ax, xlim, ylim, npix=(500,500)):
        mask = self.root.get_map(xlim=xlim, ylim=ylim, npix=npix)
        cmap = mpl.colormaps['binary']
        extent = np.append(xlim, ylim)
        ax.imshow(mask, cmap=cmap, alpha=0.5, vmin=0, vmax=3,
                  origin='lower', extent=extent)

    def plot_cell(self, ax, **kwargs):
        if 'limits' in kwargs:
            xlim, ylim = kwargs['limits']
        else:
            xlim, ylim = self.get_limits()
        if 'npix' in kwargs:
            npix = kwargs['limits']
        else:
            npix = (500,500)
        self.plot_cell_map(ax, xlim, ylim, npix=npix)

    def plot_trace(self, ax):
        """Plot incident and output rays."""
        root = self.root
        lines = {'in':'b-','out':'r-'}
        
        # select subset of events for plotting 
        if self.traces is None:
            r = self.r
            vi = self.vi
            vf = self.vf
            is_in = self.root.is_in
        else:
            r = self.r[:,self.traces]
            vi = self.vi[:,self.traces]
            vf = self.vf[:,self.traces]
            is_in = self.root.is_in[self.traces]
        
        
        # collect all incident rays
        for path in ['in','out']:
            rays = []
            if self.traces is None:
                _min,_mout = root._isect[path]['mask']
                _tin,_tout = root._isect[path]['times']
            else:
                _min = root._isect[path]['mask'][0][:,self.traces]
                _tin = root._isect[path]['times'][0][:,self.traces]
                _mout = root._isect[path]['mask'][1][:,self.traces]
                _tout = root._isect[path]['times'][1][:,self.traces]               
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
                                rin = (r[:,i] + vi[:,i]*_tin[j,i])[0::2]
                                if q_out[j,i]:
                                    rout = (r[:,i] + vi[:,i]*_tout[j,i])[0::2]
                                else:
                                    rout = r[:,i][0::2]
                                rays.append([rin,rout])
                        else:
                            if q_out[j,i]:
                                if q_in[j,i]:
                                    rin = (r[:,i] + vf[:,i]*_tin[j,i])[0::2]
                                else:
                                    rin = r[:,i][0::2]
                                rout = (r[:,i] + vf[:,i]*_tout[j,i])[0::2]
                                rays.append([rin,rout])
            if len(rays)>0:
                rays = asnumpy(cp.asarray(rays))
                cross_in = rays[:,0,:]
                cross_out = rays[:,1,:]
                ax.plot(cross_in[:,0], cross_in[:,1], 'b.', markersize=2)
                ax.plot(cross_out[:,0], cross_out[:,1], 'r.', markersize=2)
                for i in range(rays.shape[0]):
                    x = rays[i,:,0]
                    y = rays[i,:,1]
                    ax.plot(x, y, lines[path], linewidth=0.5)


    def plot(self, title='', trace=True, limits=None):
        if limits:
            xlim, ylim = limits
        else:
            xlim, ylim = self.get_limits()
        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.plot_cell(ax, limits=[xlim, ylim])
        if self.traces is None:
            r = self.r
            is_in = self.root.is_in
        else:
            r = self.r[:,self.traces]
            is_in = self.root.is_in[self.traces]
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
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        ax.grid()
        if title:
            plt.title(title, fontsize = 8)
        plt.show() 
    
    def run(self, tests=[]):
        for i in range(len(tests)):
            print('Test {} '.format(i+1), end='')
            try:
                tests[i]()
                print('passed')
            except Exception as e:
                print('failed')
                raise(e)
        
    
class TestSurface(TracingTest):
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
        self.options = {'composite':False, 'plot':True, 'nev':20}

    



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
        super().set_events(r, vi, vf)       
    
    def plot_surface(self, ax, surf, color):  
        # assume Cylinder || y in (x,z) plane
        o = asnumpy(surf.trg.orig[0::2,0])
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
        
    def plot_cell(self, ax, **kwargs):
        """Plot cell as a composite of surface groups, or a pixel map."""
        if self.options['composite']:
            gray = (0.5, 0.5, 0.5, 0.15)
            self.plot_group(ax, self.root, gray)
        else:
            super().plot_cell(ax, **kwargs)

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

    def get_title(self, op, inner):
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
        self.set_events(n=self.options['nev'])
        #self.plot(trace=False)
        self.root.evaluate()
        if self.options['plot']:
            self.plot(title='test moon plot')

    def test2(self):
        op=['or','or','or','and','and']; inner=[1,1,1,-1,-1]
        self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
        self.set_events(n=self.options['nev'])
        self.root.evaluate()
        if self.options['plot']:
            self.plot(title=self.get_title(op, inner))

    def test3(self):
        op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, 1]
        self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
        self.set_events(n=self.options['nev'])
        self.root.evaluate()
        if self.options['plot']:
            self.plot(title=self.get_title(op, inner))

    def test4(self):
        op=['or','or','or', 'or','and']; inner=[-1,-1,-1,-1, -1]
        self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
        self.set_events(n=self.options['nev'])
        self.root.evaluate()
        if self.options['plot']:
            self.plot(title=self.get_title(op, inner))

    def test5(self):
        op=['and','or','and', 'or','and']; inner=[-1,-1,-1,-1, -1]
        self.root = self.group_hexagon(R1=1.2, R2=2, op=op, inner=inner)      
        self.set_events(n=self.options['nev'])
        self.root.evaluate()
        if self.options['plot']:
            self.plot(title=self.get_title(op, inner))


    def test0(self):
        self.root = self.group_moon()
        xlim, ylim = self.get_limits()
        npix=(500, 500)
        mask = self.root.get_map(xlim=xlim, ylim=ylim, npix=npix)
        
        fig, ax = plt.subplots(figsize=(5,5))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        cmap = mpl.colormaps['binary']
        ax.imshow(mask, cmap=cmap, alpha=0.5, vmin=0, vmax=2,
                  origin='lower', extent=xlim+ylim)
        ax.grid()
        plt.title('test surface plot', fontsize = 8)
        plt.show() 
    
    def run(self, tests='all', n=20, plot=True, composite=False):
        """Test tracing through groups of surfaces.
        
        Parameters
        ----------
        tests : str or list
            Which test to run, either 'all' or a list of indices (0..5)
        n : int
            Number of events to show
        plot : bool
            Show plot with traces and intersections.
        composite : bool
            Show the cell as a composition of sub-groups with different colors.
            Otherwise, the cell is shown as a single gray area.        
        """
        self.options['composite'] = composite
        self.options['plot'] = plot
        self.options['nev'] = n
        testfnc = [self.test1, self.test2, self.test3, 
                   self.test4, self.test5]
        if isinstance(tests,list):
            torun = testfnc[tests]
        else:
            torun = testfnc
        super().run(tests=torun)


class TestConv(TracingTest):
    def __init__(self):
        super().__init__()
        self.cell = None
        
    #override 
    def get_limits(self):
        r = asnumpy(self.integ._ev['r'])
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
    
    def initiate(self):
        # load sampling events
        evdata = load_data('events_S_1mm.dat',kind='data')
        my_data = {'data':evdata,'columns':[1, 4, 7, 10, 11]}
        self.sampling = StrainSampling(my_data)
        self.sampling.print_properties()
        # define material with beam attenuation
        att = load_data('mu_Cu.dat',kind='tables')
        self.mat = Material()
        self.mat.set_extinction(Extinction.TABLE, table=att)
   
    
    def _get_traces(self, n=20, steps=None):
        idx = np.linspace(0, n, num=n, endpoint=False, dtype=int)
        if steps is None:
            nsteps = 21
            npos = 5
            steps = [(1+int(nsteps/npos))*i for i in range(npos)]
        ns = len(steps)
        res = np.zeros(n*ns, dtype=int)
        for i in range(ns):
            istp = steps[i]
            res[i*n:(i+1)*n] = idx + istp*self.integ.nev
        return res                
    
    def plot_conv(self, conv):
        if self.scan.shape  == ScanGeometry.CIRCULAR:
            ftit = 'Circular scan'
        else:
            ftit = 'Linear scan'
        xtit = 'scan position, {}'.format(self.scan.units)
        ytit = 'centre of gravity, {}'.format(self.scan.units)
        
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6))
        xmin = min(conv['x'])
        dx = max(conv['x'])-min(conv['x'])
        xlim = {'left':xmin-0.05*dx, 'right':xmin+1.05*dx}
        
        # counts
        ax=axs[0,0]
        ax.errorbar(conv['x'], conv['cnts'], yerr=conv['err_cnts'], fmt='ko-', 
                    markersize=3)
        ax.set_xlim(**xlim)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit)
        ax.set_ylabel('counts, rel. units')
        ax.grid()
        
        # pseudo-strain
        ax=axs[0,1]
        ax.errorbar(conv['x'], 1e6*conv['eps'], yerr=1e6*conv['err_eps'], 
                    fmt='ko-', markersize=3)
        ax.set_xlim(**xlim)
        #ax.set_ylim(bottom=0)
        ax.set_xlabel(xtit)
        ax.set_ylabel('pseudo-strain, 1e-6')
        ax.grid()
        
        # scan position shift and width
        ax=axs[1,0]        
        y = conv['pos'] - conv['x']        
        ax.errorbar(conv['x'], y, yerr=conv['err_pos'], fmt='bo-',  
                    markersize=3, label='shift')
        ax.errorbar(conv['x'], conv['width'], fmt='rx-',  
                    markersize=3, label='width')
        ax.set_xlim(**xlim)
        ax.legend()
        ax.set_xlabel(xtit)
        ax.set_ylabel(ytit)
        ax.grid()
        
        # cog
        ax=axs[1,1]
        ax.errorbar(conv['x'], conv['cog'][0], fmt='ro-', markersize=3, 
                    label='x')
        ax.errorbar(conv['x'], conv['cog'][1], fmt='go-', markersize=3, 
                    label='y')
        ax.errorbar(conv['x'], conv['cog'][2], fmt='bo-', markersize=3, 
                    label='z')
        if self.scan.shape == ScanGeometry.CIRCULAR:
            # add radial component for circular scan
            scan_pos = to_numpy(self.scan.to_scan_coord(conv['cog']))
            ax.errorbar(conv['x'], scan_pos['r'], fmt='ko-', markersize=3, 
                        label='r')
        ax.legend()
        ax.set_xlabel(xtit)
        ax.set_ylabel('centre of gravity, mm')
        ax.grid()
        # plot
        fig.suptitle(ftit)
        fig.tight_layout()
        plt.show()

    def test1(self, nev=30000):
        "Test linear scan."
        # define cell
        sf1 = {'R':10, 'inner':1, 'orig':[0.0, 0.0, 0.0]}           
        sf2 = {'R':10, 'orig':[10, 0.0, 0.0]}
        self.root = self.create_group(surfaces=[sf1, sf2], op='and', orig=[-10,0,0])
        self.cell = Cell(self.root, self.mat)
        self.cell.reference = 0
        # define scan
        nsteps = 21
        steps = np.linspace(-1, 11, num=nsteps)
        self.scan = ScanGeometry.linear(steps, axis=[1,0,0])
        
        # define integral and events x scan steps 
        self.integ = GaugeIntegral(self.cell, self.sampling, self.scan)
        self.integ.initiate(nev=nev)
        # define steps for event maps
        nshow = 20
        npos = 5
        steps = [(1+int(nsteps/npos))*i for i in range(npos)]
        self.traces = self._get_traces(n=nshow, steps=steps)
        # set events to the root cell
        self.set_events(self.integ._ev['r'], 
                        self.integ._ev['ki'], 
                        self.integ._ev['kf'])
        # plot event map
        self.plot()
        # make convolution integral 
        conv = self.integ.conv_func(None)
        # plot results
        self.plot_conv(to_numpy(conv))
        #for k in conv:
        #    print('{}:\n{}'.format(k,conv[k]))

    def test2(self, nev=30000, comp='rad'):
        """Test circular scan for rad or hoop component."""
        ang_comp = {'rad':0, 'hoop':90}
        # define cell
        # angle between rods
        dom = 60
        co = np.cos(dom*np.pi/180)
        si = np.sin(dom*np.pi/180)
        # radial position of rods
        a = 11
        # sample rotation
        om = 45 + ang_comp[comp]*90 
        com = np.cos(om*np.pi/180)
        som = np.sin(om*np.pi/180)
        sf1 = {'R':16, 'inner':-1, 'orig':[0.0, 0.0, 0.0]}
        sf2 = {'R':3, 'inner':1, 'orig':[a, 0.0, 0.0]}
        sf3 = {'R':3, 'inner':1, 'orig':[a*co, 0.0, a*si],'angles':[0,dom,0]}
        sf4 = {'R':3, 'inner':1, 'orig':[a*co, 0.0, -a*si],'angles':[0,-dom,0]}
        self.root = self.create_group(surfaces=[sf1, sf2, sf3, sf4], op='and', 
                                      orig=[-a,0,0], 
                                      angles=[0,om,0], 
                                      rotctr=[a,0,0])
        self.cell = Cell(self.root, self.mat)
        self.cell.reference = 1
        # define scan
        nsteps = 51
        steps = np.linspace(-90, 90, num=nsteps)
        origin_rot = [-a*com,0,a*som]
        origin = [0,0,0]
        if comp=='rad':
            origin = [-1.2,0,1.2]
        else:
            origin = [1.0,0,-1.0]
        self.scan = ScanGeometry.circular(steps, origin=origin, 
                                          axis=[0,1,0], 
                                          origin_rot=origin_rot, 
                                          step_size=1.0)
        
        # define integral and events x scan steps 
        self.integ = GaugeIntegral(self.cell, self.sampling, self.scan)
        self.integ.initiate(nev=nev)
        # define steps for event maps
        nshow = 20
        #steps = [(1+int(nsteps/npos))*i for i in range(npos)]
        steps = [0, int(0.25*nsteps), int(0.5*nsteps), int(0.75*nsteps), nsteps-1]
        self.traces = self._get_traces(n=nshow, steps=steps)
        # set events to the root cell
        self.set_events(self.integ._ev['r'], 
                        self.integ._ev['ki'], 
                        self.integ._ev['kf'])
        # plot event map
        o = asnumpy(self.root.trg.orig)[:,0]
        limits = [[o[0]-20, o[0]+20], [o[2]-20, o[2]+20]]
        self.plot(limits=limits)
        # make convolution integral 
        conv = self.integ.conv_func(None)
        # plot results
        self.plot_conv(to_numpy(conv))
        #for k in conv:
        #    print('{}:\n{}'.format(k,conv[k]))


    def init_bench(self):
        # define cell
        sf1 = {'R':10, 'inner':1, 'orig':[0.0, 0.0, 0.0]}           
        sf2 = {'R':10, 'orig':[10, 0.0, 0.0]}
        self.root = self.create_group(surfaces=[sf1, sf2], op='and', orig=[-10,0,0])
        self.cell = Cell(self.root, self.mat)
        self.cell.reference = 0
        # define scan
        nsteps = 25
        steps = np.linspace(-1, 11, num=nsteps)
        self.scan = ScanGeometry.linear(steps, axis=[1,0,0])
        
        # define integral and events x scan steps 
        self.integ = GaugeIntegral(self.cell, self.sampling, self.scan)
        self.integ.initiate(nev=40000)
        # set events to the root cell
        self.set_events(self.integ._ev['r'], 
                        self.integ._ev['ki'], 
                        self.integ._ev['kf'])
        
        
    def bench_func(self):
        #f = np.random.random()
        self.integ.conv_func(None)

    def run(self):
        self.initiate()
        super().run(tests=[self.test1, self.test2])
        #super().run(tests=[self.test2])
        
    def bench(self, n_repeat=10):
        if has_gpu():
            from cupyx.profiler import benchmark
            #self.initiate()
            self.init_bench()
            res =  benchmark(self.bench_func, (), n_repeat=n_repeat)
            print(res)
            time = cp.average(res.gpu_times) + cp.average(res.cpu_times)
        else:
            import timeit
            glob = globals()
            glob.update({'self':self})
            res = timeit.timeit('self.bench_func()',
                                globals=glob, 
                                setup='self.init_bench()',
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

if __name__ == "__main__":
    print('running main')
    if 'primitives' in _run_tests:
        test = TestSurface()
        test.run(n=100, plot=True, composite=False)
    if 'integrals' in _run_tests:
        test = TestConv()
        test.run()
    if 'bench' in _run_tests:
        bench_convolution()
        


