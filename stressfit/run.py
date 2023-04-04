# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""


Created on Thu Mar 16 11:14:38 2023
@author: saroun
"""
import stressfit.dataio as dataio
from stressfit.ui.config import uiconfig
import stressfit.commands as comm
import stressfit.graphs as gr
import stressfit.shapes as shapes
import numpy as np

def _init(obj):
    """Initialize StressFit run configuration.
    
    Sets attributes with references to the objects handling various tasks
    in StressFit application:
        
    - logger : Handles loggin tasks. 
    - workspace : handles workspace setting.
    - config : Handles user input data.
    - input : Permits to load/save complete user input in json file.
        
    When creating these handlers, initialization tasks are performed such as 
    the default workspace setting, loading of default input data, etc.
    
    Parameters
    ----------
    obj : obj
        An object for which this function sets the attributes.
    """
    obj.logger = dataio.logger()
    # Set workspace data and create workspace directories if not yet done.
    obj.workspace  = dataio.workspace()
    # Create and load program confguration. 
    obj.config  = uiconfig()
    # Program input handler
    obj.input = _InputFileHnd()

class _InputFileHnd():
    """Handles save/load of program input files."""
    
    def __init__(self):
        self._last_input_file = 'input.json'
        self._conf = uiconfig()
        self._log = dataio.logger()
    
    def import_from(self, src, json=True):
        self._conf.import_from(src, json=json, reload=True)
        if not self._conf.is_ready():
            self._log.error('Program input is not complete.') 
    
    def load(self, filename='', path='', verbose=True):
        """Load program input file (json format)."""
        if not filename:
            filename = self._last_input_file
        try:
            lines = dataio.load_text(filename, kind='work', path=path,
                                     verbose=verbose)
            self._last_input_file = filename
            # import into global input data
            self._conf.import_from(lines, reload=True)
            if not self._conf.is_ready():
                self._log.error('Program input is not complete.')
                return            
        except Exception:
            msg = 'Cannot load program configuration {}.'
            self._log.exception(msg.format(filename))

    def save(self, filename='', path=''):
        """Save program input file (json format)."""
        if not filename:
            filename = self._last_input_file
        fn = dataio.get_input_file(filename, kind='work', path=path)
        try:
            out = self._conf.export_to()
            dataio.save_text(out, fn)
            self._last_input_file = filename
        except Exception:
            msg = 'Cannot save program configuration in {}'
            self._exception(msg.format(fn))


class _Fitting():
    
    def __init__(self, name='fit_imodel'):
        if not name in ['fit_imodel', 'fit_emodel']:
            raise Exception('Unknown fitting task: {}'.format(name))
        self.logger = dataio.logger()
        self.conf =  uiconfig()
        self.name = name
        if name == 'fit_imodel':
            self._cid = 'imodel'
        else:
            self._cid = 'emodel'

    def _create_fitobj(self, scan, model, scale, par):
        if self._cid == 'imodel':
            fitobj = comm.define_ifit(scan, model['dist'], par['nrec'], 
                                    ndim=par['npts'])        
            fitobj.defScaling(scale['values'], scale['fit'], 
                            minval=[0., 0., -np.inf])
            fitobj.setInterpModel(model['interp'])
        else:
            fitobj = comm.define_sfit(scan, model['dist'], par['nrec'], 
                        ndim=par['npts'], 
                        z0=scale['z0'],
                        eps0=scale['eps0'])  
        fitobj.setInterpModel(model['interp'])
        return fitobj        

    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        pfx = ''
        par = self.conf.get_config(self.name)
        dname = par['data']
        if pfx:
            base = '{}_{}_{}'.format(self._cid, pfx, dname)
        else:
            base = '{}_{}'.format(self._cid, dname)
        if ext:
            fname = base + ext
        else:
            fname = base
        return fname
    
    def _prepare(self):
        """Make all settings necessary to run commands.
        
        Returns
        -------
        MCCFit object
        
        """
        # reload all input data if needed
        self.conf.reload_all()
        # get command parameters
        par = self.conf.get_config(self.name)
        # selected scan data 
        scan = self.conf.get_scan(par['data'])
        # selected model
        model = self.conf.get_item(self._cid, item=par['model'])
        # scale parameters
        scale = model['scale']
        # set scan parameters and sampling
        comm.set_scan(scan)        
        # set attenuation
        att = self.conf.attenuation
        comm.set_attenuation(att)        
        # create fit object
        fitobj = self._create_fitobj(scan, model, scale, par)
        return fitobj

    def _on_replot(self):
        try:
            fitobj = self._prepare()        
            # create smeared curve: run without fitting, maxiter=0
            comm.run_fit(fitobj, maxiter=0, outname='')
        except Exception as e:
            self.logger.exception(str(e))
        
    def _on_guess(self):
        try:
            fitobj = self._prepare()    
            par = self.conf.get_config(self.name)
            comm.run_fit_guess(fitobj,
                               maxiter=par['fit']['maxiter'], 
                               ar=par['fit']['ar'],
                               outname='')
        except Exception as e:
            self.logger.exception(str(e))
            
    def _on_fit(self, save=True):
        try:
            fitobj = self._prepare()  
            par = self.conf.get_config(self.name)             
            comm.run_fit(fitobj, 
                         maxiter=par['fit']['maxiter'], 
                         loops=par['fit']['loops'],
                         ar=par['fit']['ar'],
                         guess=par['fit']['guess'],
                         bootstrap=par['fit']['loops']>2,
                         outname=None)
            if save:
                fname = self._get_output_filename()
            else:
                fname = ''
            comm.report_fit(fitobj, fname)
        except Exception as e:
            self.logger.exception(str(e))            

    def _on_reg(self, save=True):
        try:
            fitobj = self._prepare()  
            par = self.conf.get_config(self.name)
            rang = par['reg']['range']
            steps = par['reg']['steps']
            dr = (rang[1]-rang[0])/(steps-1)
            ar = steps*[0]
            for i in range(steps):
                ar[i] = rang[0] + i*dr
            if save:
                fname = self._get_output_filename()
            else:
                fname = ''
            comm.run_fit_reg(fitobj, 
                             maxiter=par['fit']['maxiter'], 
                             ar=ar, 
                             guess=par['fit']['guess'],
                             outname=fname,
                             )
        except Exception as e:
            self.logger.exception(str(e)) 


class _Runner():
    TASKS=['scene',
           'strain',
           'gauge', 
           'data', 
           'intensity_fit',
           'strain_fit',
           'strain_reg']
    
    def __init__(self, interactive=True, work=None, config=None):
        import matplotlib
        _init(self)
       # self.logger = dataio.logger()
        # Set workspace data and create workspace directories if not yet done.
       # self.workspace = dataio.workspace()
        # Create and load program confguration. 
        # Reload all input data defined by this configuration
       # self.config = uiconfig()
        # Program input handler
       # self.input = _InputFileHnd()
        self.shapes = dataio.load_config('shapes')
        
        isi = matplotlib.is_interactive()
        # switch to non-interactive mode
        if isi and not interactive:
            matplotlib.use('agg')
        # self._plt_backend = matplotlib.get_backend()
        if work and isinstance(work, str):
            self.workspace.root = work
            self.workspace.info()
        if isinstance(config, str):
            if config.startswith('{'):
                self.input.import_from(config)
            else:
                self.load(file=config)
        elif isinstance(config, dict):
            self.input.import_from(config, json=False)
    
    def load(self, file='', path=''):
        if file:
            self.input.load(filename=file, path=path)

    def _task_to_cid(self, task):
        cid = None
        what = None
        if task in _Runner.TASKS:
            cid = task
            if task=='strain':
                cid = 'resolution'
                what = ['strain']
            elif task=='gauge':
                cid = 'resolution'
                what = ['resolution']
            elif task=='intensity_fit':
                cid = 'fit_imodel'
                what = ['fit']
            elif task=='strain_fit':
                cid = 'fit_emodel'
                what = ['fit']
            elif task=='strain_reg':
                cid = 'fit_emodel'
                what = ['reg']
                
        elif task in self.config.config_keys():
            cid = task
            if cid=='resolution':
                what = ['strain',  'resolution']
            elif cid=='fit_imodel':
                what = ['plot', 'fit']
            elif task=='fit_emodel':
                what = ['plot', 'fit']
            elif task=='fit_emodel':
                what = ['reg']
        return cid, what
        
    
    def _set_shape(self):
        try:
            par = self.config.get_item('shape')
            if par['shape'] == 'File':
                fn = par['param']['filename']
                fullname = dataio.get_input_file(fn['file'], path=fn['path'])
                obj = shapes.create(par['shape'], filename=fullname)
            else:
                obj = shapes.create(par['shape'], **par['param'])        
            comm.set_shape(obj)
        except:
            self.logger.exception('Cannot create sample shape object.')
    
    def define_shape(self, name, **kwargs):
        """Define new sample shape.
        
        Run :func:`stressfit.shapes.help` for available shape names and 
        parameters.
        """
        if not name in self.shapes:
             self.logger.error('Unknown shape name: {}'.format(name))
             return
        par = self.config.get_item('shape')
        new_param = self.shapes[name]['param']
        par['param'].clear()
        for key in new_param:
            par['param'][key] = new_param[key]['value']
        
        par = self.config.get_item('shape')
        if name != par['shape']:
            new_param = self.shapes[name]['param']
            par['shape'] = name
            par['param'].clear()
            for key in new_param:
                par['param'][key] = new_param[key]['value']
        if len(kwargs)>0:
            self.config.set_item('shape', {'shape': name, 'param': kwargs})
    
    def info_shape(self):
        """Print information on current sample shape definition."""
        data = self.config.get_item('shape')
        self.logger.info('\nSample shape: {}'.format(data['shape'])) 
        self.logger.info('\tparam={}'.format(data['param']))
    
    def info(self):
        """Print information on current run configuration data."""
        ver = dataio.__package_info__['version']
        self.logger.info('\nStressFit ver={}'.format(ver))
        fmt = 20*'-'
        self.logger.info(fmt)
        self.list_tasks()
        self.list_shapes()
        self.info_shape()
        self.list_geom()
        self.list_sampling()
        self.list_data()
    
    def list_data(self):
        """List information on data sets defined in the run configuraiton."""
        self.logger.info('\nData:')  
        data = self.config.get_item('data')
        fmt = '\t{}:\n{}'
        rec_keys = ['strain','intensity','geometry','sampling']
        fmt_rec = '\t\tstrain={}, intensity={}, geometry={}, sampling={}'
        for key in data:
            rec = [data[key][k] for k in rec_keys]
            rec_str = fmt_rec.format(*rec)
            s = fmt.format(key,rec_str)
            self.logger.info(s)

    def list_geom(self):
        """List information on geometries defined in the run configuraiton."""
        self.logger.info('\nGeometries:') 
        data = self.config.get_item('geometry')
        fmt = '\t{}:\n{}'
        rec_keys = ['angles','rotctr','scandir','scanorig']
        fmt_rec = '\t\tangles={}, rotctr={}, scandir={}, scanorig={}'
        for key in data:
            rec = [data[key][k] for k in rec_keys]
            rec_str = fmt_rec.format(*rec)
            s = fmt.format(key,rec_str)
            self.logger.info(s)

    def list_sampling(self):
        """List information on sampling data defined in the run configuraiton."""
        self.logger.info('\nSampling:') 
        data = self.config.get_item('sampling')
        fmt = '\t{}:\n{}'
        fmt_rec = '\t\tfile={}, nrec={}, loaded={}'
        for key in data:
            fd = data[key]['fdata']
            rec = [data[key]['file']['file'], data[key]['nrec'], fd is not None]
            rec_str = fmt_rec.format(*rec)
            s = fmt.format(key,rec_str)
            self.logger.info(s)
            
    def list_tasks(self):
        """Print a list of defined tasks.
        
        Tasks to be used as an argument to the methods 
        :meth:`run`, :meth:`set_par` and :meth:`get_par`.
        """
        self.logger.info('\nTasks:')
        fmt = '\t{}'
        for key in _Runner.TASKS:
            self.logger.info(fmt.format(key))
        
    def list_shapes(self):
        """Print a list of available shapes.
        
        For complete information, run :meth:`stressfit.shapes.help`.
        """
        self.logger.info('\nShapes:')
        fmt = '\t' + (len(self.shapes)-1)*'{}, ' + '{}'
        k = list(self.shapes.keys())
        out = fmt.format(*k)
        self.logger.info(out)
            
    def list_param(self, task):
        """List parameters for given task.
        
        Use :meth:`get_tasks` for a list of defined tasks.
        
        """
        cid, what = self._task_to_cid(task)
        par = self.config.get_config(cid)
        self.logger.info('\nParameters for {}:'.format(task))
        fmt = '\t{} = {}'
        for key in par:
            self.logger.info(fmt.format(key, par[key]))

    def set_param(self, task, **kwargs):
        """Set parameters for given task.
        
        Use :meth:`get_tasks` for a list of defined tasks.
        
        Parameters
        ----------
        task : str
            Task name.
        kwargs : dict
            Input data provided as dict. Use :meth:`get_par` to get
            the current data for given task.
        
        """
        try:
            if len(kwargs) > 0:
                cid, what = self._task_to_cid(task)
                self.config.set_config(cid, kwargs)
        except:
            msg = 'Wrong input parameters for {}: {}'
            self.logger.exception(msg.format(task, kwargs)) 
         
    
    def run(self, task, save=True, **kwargs):
        if not task in _Runner.TASKS:
            self.logger.error('Undefined task: {}'.format(task))
            return
        cid, what = self._task_to_cid(task)
        if 'what' in kwargs:
            what = kwargs['what']
            del kwargs['what']
        self._set_shape()
        self.logger.info('Running task: {}'.format(task))
        if task=='scene':
            self._plot_scene(save=save, **kwargs)
        elif task=='strain':
            self._plot_strain(save=save, what=['strain'], **kwargs)
        elif task=='gauge':
            self._plot_strain(save=save, what=['resolution'], **kwargs)
        elif task=='data':
            self._plot_data(save=save, **kwargs)
        elif task=='intensity_fit':
            self._fit_intensity(save=save, what=what, **kwargs)
        elif task=='strain_fit':
            self._fit_strain(save=save, what=what, **kwargs)
        elif task=='strain_reg':
            self._fit_strain(save=save, what=what, **kwargs)
         
    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        par = self.config.get_config('scene')
        pfx = ''
        ori = par['geometry']
        spl = par['sampling'] 
        if pfx:
            base = '{}_{}_{}'.format(pfx, spl, ori)
        else:
            base = '{}_{}'.format(spl, ori)
        if ext:
            fname = base + ext
        else:
            fname = base
        return fname        
    
    def _get_output_filename_data(self, name='data_name', ext=''):
        """Generate base output filename."""
        pfx = ''    
        if pfx:
            base = '{}_{}'.format(pfx, name)
        else:
            base = '{}'.format(name)
        if ext:
            fname = base + '_{}'.format(ext)
        else:
            fname = base
        return fname 
    
    
    def _plot_scene(self, save=True, **kwargs):
        # get command parameters
        par = self.config.get_config('scene')
        self.config.reload('sampling', item=par['sampling'])
        if not self.config.is_ready():
            self.logger.error('Input data not ready.')
            return 
        self.set_param('scene', **kwargs)
        # set selected sampling and geometry  
        g = self.config.get_item('geometry',item=par['geometry'])
        s = self.config.get_item('sampling',item=par['sampling'])['fdata']
        comm.set_geometry(g)
        comm.set_sampling(s)
        # do plot
        if save:
            fname = self._get_output_filename(ext='.png')
        else:
            fname = ''
        comm.plot_scene(par['nrec'], 
                        filename=fname, 
                        rang=2*[par['rang']],
                        proj=par['proj'])
    
    def _plot_strain(self, save=True, what=['strain', 'resolution'], **kwargs):       
        # get command parameters
        par = self.config.get_config('resolution')
        self.config.reload('sampling', item=par['sampling'])
        self.config.reload('attenuation')
        if not self.config.is_ready():
            self.logger.error('Input data not ready.')
            return
        self.set_param('resolution', **kwargs)
        # set selected sampling and geometry
        g = self.config.get_item('geometry',item=par['geometry'])
        s = self.config.get_item('sampling',item=par['sampling'])['fdata']
        comm.set_geometry(g)
        comm.set_sampling(s)
        # set attenuation
        att = self.config.attenuation
        comm.set_attenuation(att)

        # do plot
        nrec = par['nrec']
        fname = self._get_output_filename()
        rang = list(par['rang'])
        nstp = par['steps']
        scan_range = rang + [nstp]
        if isinstance(what, str):
            what = [what]
        try:
            for key in what:
                if key=='strain':
                    comm.report_pseudo_strains(scan_range, fname, 
                                               nev=nrec,
                                               intensity=True,
                                               inline=True,
                                               plot=True, 
                                               save=save)
                if key=='resolution':    
                    comm.report_resolution(scan_range, fname, 
                                           nev=nrec,
                                           cog=True,
                                           inline=True,
                                           plot=True, 
                                           save=save)
        except Exception as e:
            self.logger.exception(str(e))

    def _plot_data(self, save=True, **kwargs):
        """Plot data with simulated pseudo-strains and pseudo-intensities."""
        # get command parameters
        par = self.config.get_config('data')
        self.config.reload_all()
        if not self.config.is_ready():
            self.logger.error('Input data not ready.')
            return
        self.set_param('data', **kwargs)
        # set attenuation
        att = self.config.attenuation
        comm.set_attenuation(att)       
        # collect simulated pseudo-strains and experimental data to plot 
        dlist = self.config.get_item('data')
        expdata = {}
        simdata = {}
        for name in dlist:
            scan = self.config.get_scan(name)
            x = scan['eps'][:,0]        
            scan_range = [min(x), max(x), 2*len(x)+1]
            fname = self._get_output_filename_data(name=scan['epsfile'])       
            nrec = par['nrec']
            save = save
            # set geometry and sampling for given scan
            comm.set_geometry(scan)
            comm.set_sampling(scan['sampling'])
            res = comm.report_pseudo_strains(scan_range, fname, 
                                             nev=nrec,
                                             intensity=True,
                                             inline=True,
                                             plot=False, 
                                             save=save)
            expdata[name] = scan
            simdata[name] = res
        # do plot
        gr.plot_comparison(simdata, expdata, 
                           title='Experimental data vs. pseudo-stran')        
        
    def _fit_intensity(self, save=True, what=['plot', 'fit'], **kwargs):
        fitter = _Fitting(name='fit_imodel')
        self.logger.info('Fitting of intensity model.')
        self.set_param('fit_imodel', **kwargs)
        if isinstance(what, str):
            what = [what]
        for key in what:
            if key=='plot':
                fitter._on_replot()
            if key=='fit':
                fitter._on_fit(save=save)
            if key=='reg':
                fitter._on_reg(save=save)
        
    def _fit_strain(self, save=True, what=['plot', 'fit'], **kwargs):
        fitter = _Fitting(name='fit_emodel')
        self.logger.info('Fitting of strain model.')
        self.set_param('fit_emodel', **kwargs)
        if isinstance(what, str):
            what = [what]
        for key in what:
            if key=='plot':
                fitter._on_replot()
            if key=='fit':
                fitter._on_fit(save=save)
            if key=='reg':
                fitter._on_reg(save=save)        

def run(work=None, config=None, tasks=['scene', 'strain', 'gauge'],
        info=False):
    """Run given tasks with stressfit.
    
    Provide tasks either as list or dict. List provides task
    names to be excuted with current input data. A dictionary would provide
    task names as keys, and task parameters as values. The values should
    be a dictionary with parameters to be changed before the task is executed.
    
    Parameters
    ----------
    work : str
        Workspace root directory.
    load : str
        Input configuration to be loaded (json format)
    tasks : dict or list
        Task names with optional parameters.
    
    """
    r = _Runner(work=work, config=config)
    
    if isinstance(tasks, list):
        if info: r.info()
        for task in tasks:
            r.run(task)
    elif isinstance(tasks, dict):
        if info: r.info()
        for task in tasks:
            param = tasks[task]
            if isinstance(param, dict):
                r.run(task, **param)
            else:
                msg = 'Task {} has wrong parameters format. Using default.'
                r.logger.warning(msg.format(task))
                r.run(task)
    else:
        r.info()
        r.logger.info('No tasks specified for running.')

#%% test

def test_1():
    """Test _Runner, basic."""
    r = _Runner()
    r.info()
    task = 'strain_fit'
    r.set_param(task, fit={'maxiter':100, 'ar':3}, nrec=3000)
    r.list_param(task)
    #r.run(task)

def test_3(work='', load='', task='', **kwargs):
    """Test _Runner, switching to test workspace and loads configuraiton.
    
    Parameters
    ----------
    work : str
        Workspace root directory.
    load : str
        Input configuration to be loaded (json format)
    task : str
        Task name (use _Runner.info() for the list of tasks)
    kwargs : dict
        Input parameters for given task.
    
    """
    r = _Runner(work=work, config=load)
    r.info()
    if task:
        r.run(task, **kwargs)
        r.list_param(task)
    #r.run(task)


def test_2(tasks=['scene', 'strain', 'gauge']):
    """Test _Runner for given execution tasks."""
    run(tasks=tasks)
    
if __name__ == "__main__":
    import sys
    keys = ['work', 'config', 'task']
    args = sys.argv[1:]
    #args = ['work=D:\Saroun\Publikace\2023\ECNS2023\stressfit - Copy', 
    #        'config=input_tube_scan_HK4.json', 'gauge']
    tasks = []
    kw = {}
    for a in args:
        if a.find('=')>0:
            cm = a.split('=')
            if cm[0] in keys:
                if cm[0] == 'task':
                    tasks.append(cm[1])
                else:
                    kw[cm[0]] = cm[1]
        elif a in _Runner.TASKS:
            tasks.append(cm[0])
    if len(tasks)>0:
        kw['tasks'] = tasks
    kw['info'] = True
    print(kw)
    run(**kw)
    
            
    
    