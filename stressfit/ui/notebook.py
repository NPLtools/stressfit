"""Classes and functions for creating notebook interface.

Created on Tue Oct  5 11:38:15 2021
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import abc
import numpy as np
import ipywidgets as ipy
import json
from IPython.display import display
from .widgets import DirInput, FileInput, ArrayInput, SelectInput
from .widgets import create_header, create_select, create_input_int 
from .widgets import create_input_float
from .widgets import SButton, SRadioButtons, SCheckbox

import stressfit.commands as comm
import stressfit.shapes as shapes
import stressfit.dataio as dataio

def _has_keys(args, keys):
    """Verify that all args are in keys."""
    return all (k in args for k in keys)

#%% Input collections

class UI_base:
    """Common parent class for all input collection."""
    
    def _err_ni(src):
        msg = "{}(): Subclass must implement abstract method".format(src)
        raise NotImplementedError(msg)
    
    _NON_IMPLEMENTED = "Subclass must implement abstract method"
    def __init__(self, name, keys):
        self._err = None # error output, can be defined by ui()
        self._msg = None # message output, can be defined by ui()
        self._name = name
        self._keys = keys # list of allowed value keys
        self._widgets = {} # dict of input widgets
        self._values = {} # dict of widget values
        self._on_change = None
        self._on_init = None
    
    def _check_keys(self,data):
        """Verify that data contains all required keys."""
        if not _has_keys(data,self._keys):
            msg = 'Required key(s) missing in the input data: {}'
            raise Exception(msg.format(self._keys))
    
    @property
    def name(self):
        """ID string for the input block."""
        return self._name
            
    ### Abstract methods to be overriden
   
    @abc.abstractmethod
    def update_widgets(self, data):
        """Update input widgets from data.

        Parameters
        ----------
        data: dict
            Input names and values 
        
        """        
        UI_base._err_ni('update_widgets')

    @abc.abstractmethod
    def update_values(self):
        """Update input data from widgets.
        
        Input data are stored in self._values as dict.
        
        """        
        UI_base._err_ni('update_values')

    ### Other methods
    def show(self, widgets=[], err=None, msg=None, **kwargs):
        """Display VBox container with all input widgets."""        
        self._err = err
        self._msg = msg
        display(ipy.VBox(widgets))

    def message(self, txt, clear=True):
        """Print message in the message output, if defined."""
        if self._msg:
            if clear:
                self._msg.clear_output()
            with self._msg:
                print(txt)
                
    def error(self, txt, clear=True):
        """Print error message in the error output, if defined."""
        if self._err:
            if clear:
                self._err.clear_output()
            with self._err:
                s = "<font color='red'>{}</font>".format(txt)
                t = ipy.HTML(value=s)
                display(t)
                
    def get_values(self):
        """Return actual input values as dict."""
        self.update_values()
        return self._values

    def set_on_change(self, on_change):
        """Set callback to be executed when the UI changes state."""
        self._on_change = on_change
        
    def _call_change(self, **kwargs):
        """Call when the UI changes state."""
        if self._on_change:
            self._on_change(self, **kwargs)
    
    def set_on_init(self, _on_init):
        """Set callback to be executed when the UI needs to initialize."""
        self._on_init = _on_init
        
    def _call_init(self, **kwargs):
        """Call when the UI changes state."""
        if self._on_init:
            self._on_init(self, **kwargs)

class UI_workspace(UI_base):
    """Input for handling workspace directories.
     
    Defines the root workspace directory and several user directories,
    which can be either absolute or relatve to the workspace root.
    The workspace is encapsulated by a dataio.Workspace object.
     
    The user directory types are:
        work:
            Root workspace directory   
        data
            Path to input data
        tables
            Path to material tables etc.
        instruments
            Path to instrumemnt definitions
        output
            Path to output directory
    
    When the work directory is changed to a new one with saved workspace 
    configuration, the other directories are updated.

    In addition, this UI enables to save, reload and reset the configuration.
    If reset is chosen, the package default is restored.
     
    Parameters
    ----------
    exclude: list
        List path types which should be excluded from the input list.

    """
    
    _inp = {'work':['Workspace', 'Root workspace directory'],
            'data':['Input', 'Path to input data'],
            'tables':['Tables', 'Path to material tables etc.'],
            'instruments':['Instruments', 'Path to instrumemnt definitions'],
            'output':['Output', 'Path to output directory']}
    
    def __init__(self, exclude = ['instruments']):
        # NOTE: the keys list must match dataio.__path_keys
        super().__init__('workspace', dataio.Workspace.types) 
        self.wk = dataio.workspace()
        self._values.update(self.wk.get_paths())
        self._out = ipy.Output(layout=ipy.Layout(width='100%'))
        for key in self._keys:
            if not key in exclude:
                [lbl, ttp] = UI_workspace._inp[key]
                inp = DirInput(name=key, 
                               path=self._values[key], 
                               label=lbl, 
                               tooltip=ttp)
                inp.observe(self._on_path_change)
                self._widgets[key] = inp
        self._widgets['work'].txt.disabled=True
        self._widgets['work'].callback = self._on_workspace_change
        
    def _on_workspace_change(self,s):
        with self._out:
            res = self.wk.change_workspace(self._values['work'])
            data = self.wk.get_paths()
            self.update_widgets(data)
            self.update_values()
            if res:
                print('Workspace configuration reloaded.')
            self._call_change({'work':data['work']})
    
    def _on_path_change(self, inp):
        key = inp['name']
        value = inp['value']
        self._values[key] = value
        self._call_change(**{key:value})
        
    def _on_apply(self, b):
        self.validate_workspace()
    
    def _on_reset(self, b):
        """Set workspace configuration to package default."""
        self.reset_workspace()
        
    def _on_save(self, b):
        """Save workspace configuration in workspace root directory."""
        self._out.clear_output()
        try:
            self.update_workspace()
            self.wk.validate_paths()
            self.wk.save()
            with self._out:
                wkp = self.wk.get_paths(keys=['work'])
                msg = 'Workspace configuration saved in {}/{}.'
                print(msg.format(wkp['work'],dataio.Workspace.cfg_name))
        except Exception as e:
            print(e)
    
    def _on_load(self, b):
        """Re-load workspace configuration from workspace root directory."""
        with self._out:
            if self.wk.load():
                data = self.wk.get_paths()
                self.update_widgets(data)
                self.update_values()
                print('Workspace configuration reloaded.')
            else:
                wkp = self.wk.get_paths(keys=['work'])
                msg = 'Workspace configuration file {} not found in {}'
                print(msg.format(dataio.Workspace.cfg_name, wkp['work']))
            
    def update_widgets(self, data):
        """Update widgets and update workspace from data.

        Parameters
        ----------
        data : dict
            Expected keys are ['work', 'data', 'tables', 'instruments', 
                               'output']
            
            All keys are optional.

        """
        for key in data:
            if key in self._widgets:
                self._widgets[key].value = data[key]

    def update_values(self):
        """Update input data from widgets.
        
        NOTE: The widgets automatically update _values so this method is
        not in fact needed. It is provided only for consistency with the other 
        UI widtget collections.
        
        """ 
        for key in self._widgets.keys():
            self._values[key] = self._widgets[key].value
    
    def show(self, err=None, msg=None, **kwargs):
        """Display widgets."""
        appl = ipy.Button(description='Apply',
                          tooltip='Update and validate workspace paths.',
                          layout=ipy.Layout(width='10%'))
        appl.on_click(self._on_apply)
        rst = ipy.Button(description='Reset',
                         tooltip='Reset workspace configuration to package default.',
                         layout=ipy.Layout(width='10%'))
        rst.on_click(self._on_reset)
        save = ipy.Button(description='Save',
                          tooltip='Save workspace configuration.',
                          layout=ipy.Layout(width='10%'))
        save.on_click(self._on_save)
        load = ipy.Button(description='Load',
                          tooltip='Load workspace configuration.',
                          layout=ipy.Layout(width='10%'))
        load.on_click(self._on_load)
        hdr = create_header('Set workspace')
        
        layout = ipy.Layout(display='flex',flex_flow='row',
                                border='none',width='100%')
        appl_box = ipy.HBox([appl, save, load, rst], layout=layout)
        
        lst = [hdr]
        for w in self._widgets:
            lst.append(self._widgets[w].ui())
        lst.append(appl_box)
        lst.append(self._out)
        super().show(widgets=lst, err=err, msg=msg)
        
    def update_workspace(self):
        """Update workspace configuration according to the widget values."""
        self._out.clear_output()
        self.wk.change_workspace(self._values['work'])
        self.wk.set_paths(**self._values)
                
    def validate_workspace(self, verbose=True):
        """Check workspace configuration (paths must exist)."""
        self._out.clear_output()
        with self._out:
            try:
                self.update_workspace()
                self.wk.validate_paths()
                if verbose: 
                    print('Workspace OK')
            except Exception as e:
                print(e)

    def reset_workspace(self):
        """Set workspace configuration to package default."""
        self._out.clear_output()
        with self._out:
            self.wk.reset_paths()
            data = self.wk.get_paths()
            self.update_widgets(data)
            self.update_values()




class UI_shape(UI_base):
    """Dialog for sample shape definition.
    
    The values are returned by the method :meth:`get_values`. It returns
    named list with the keys:
        select:
            A key for the selected sample shape.
        param:
            A named list of the selected shape parameters.
    
    Use the method :meth:`create_shape` to create an instance of selected 
    shape with given parameters.
    
    The shapes are defined in a configuration file in package repository.
    It is loaded by calling `dataio.load_config('shapes')`.
    
    Parameters
    ----------
    select: str
        Initially selected shape, e.g. stressfit.shapes.Plate.
    
    """
    
    def __init__(self, select=shapes.Plate):
        super().__init__('shape', ['shape', 'param']) 
        self.select = select
        self.shapes = dataio.load_config('shapes')
        self._values['shape'] = self.select
        self._values['param'] = self.shapes[self.select]['param']
        self.options = []
        for key in self.shapes:
            self.options.append((self.shapes[key]['name'],key))
        
        wdg = create_select(name='shape', options=self.options, 
                            label='', value=self.select,
                            width_label=0, width_drop=150)
        self._widgets['shape'] = wdg
        self._reset_widgets()
        self._widgets['shape'].observe(self._on_change_shape, type='change')
        self._out = ipy.Output(layout=ipy.Layout(border='1px solid'))

    def _reset_widgets(self):
        """Re-fill self._widgets with input widgets for selected shape."""
        param = self.shapes[self.select]['param']
        self._widgets['shape'].value = self.select
        items = {}
        for key in param:
            p = param[key]
            value =p['value'] 
            uni = p['unit'] 
            desc = p['label']
            hint = p['hint']
            if uni:
                hint = ' [{}]'.format(uni) + ' ' + hint
            if isinstance(value, float) or isinstance(value, int):
                wgt = ArrayInput(name=key, value=value, label=desc, hint=hint)
            elif isinstance(value, list):
                if isinstance(value[0],float):
                    wgt = ArrayInput(name=key, value=value, label=desc, hint=hint)
                else:
                    wgt = SelectInput(name=key, options=value, label=desc, hint=hint,
                                      width_drop=80)
            else:
                raise Exception('Unknown parameter type: {}'.format(value))
            wgt.observe(self._on_change_param)
            items[key] = wgt
        
        self._widgets['param'] = items
        self.update_values()

    def _on_change_param(self, inp):
        key = inp['name']
        if key in self._values['param']:
            self._values['param'][key] = inp['value']
            self._call_change(**{key:inp['value']})
        else:
            raise Exception('Unknown input key: {}'.format(key))

    def _on_change_shape(self, change):
        """Process change of selected shape.
        
        Create new list of widgets and display it in the output area.
        
        """
        layout = ipy.Layout(width='100%')
        if (change['name']=='value'):
            self.select = change['new']
            self._reset_widgets()
            wdg = []
            for w in self._widgets['param'].values():
                wdg.append(w.ui())
            box = ipy.VBox(wdg, layout=layout)
            self._out.clear_output()
            with self._out:
                display(box)
            self._call_change(**{'shape':self.select})
 
    def update_widgets(self, data):
        """Update parameter input widgets from data.

        Parameters
        ----------
        data: dict
            Input names and values 
        
        """
        self._check_keys(data)
        if data['shape'] != self.select:
            self._on_change_shape({'name':'value', 'new':data['select']})
        param = data['param']
        items = self._widgets['param']
        for key in param:
            item = items[key]
            if isinstance(item, SelectInput):
                items[key].set_index(param[key].index)
            else:
                items[key].value=param[key].value           
        
    def update_values(self):
        """Update input data from widgets.
        
        Input data are stored in self._values as dict.
        
        """
        self._values.clear()
        self._values['shape']=self.select
        param = {}
        items = self._widgets['param']
        for key in items:
            item = items[key]
            if isinstance(item, SelectInput):
                value = item.get_index()
            else:
                value = item.value
            param[key] = value
        self._values['param'] = param
    
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        lbl = create_header('Sample shape')
        box = [lbl, self._widgets['shape'], self._out]
        super().show(widgets=box, err=err, msg=msg)
        self._on_change_shape({'name':'value', 'new':self.select})
    
    def create_shape(self):
        """Create selected shape instance from stressfit.shapes."""
        self.update_values()
        comp = shapes.create(self._values['shape'],**self._values['param'])
        return comp

        
class UI_geometry(UI_base):
    """Dialog for sample geometry.
    
    This UI allows to set the sample and scan geometry and mainain
    a named list of such oreiantations.
    
    Each geometry is defined by four vectors:    
    
    angles : array_like(3)
        Euler YZY angles of sample geometry
    rotctr : array_like(3)
        Rotation centre in local sample coordinates [mm]
    scandir : array_like(3)
        Scan direction in local sample coordinates
    scanorig : array_like(3)
        Scan origin (position corresponding to zero scan position)
        
    The values are returned by the method :meth:`get_values`. It returns
    named list with the keys:
        input:
            Values of the above oerientation vectors from the input widgets
        list:
            A named list of oreintations added by the user. 
               
    """
    
    def __init__(self):
        super().__init__('geometry', 
                         ['angles', 'scandir', 'rotctr', 'scanorig'])   
        # widgets for current dataset
        
        inp = {}
        inp['angles'] = {'value':[135,0,0], 
                         'hint':'Sample rotation: YZY Euler angles [deg]'}
        inp['scandir'] = {'value':[0,0,1], 
                         'hint':'Scan direction vector'}
        inp['rotctr'] = {'value':[0,0,0], 
                         'hint':'Rotation centre [mm]'}       
        inp['scanorig'] = {'value':[0,0,0], 
                         'hint':'Scan origin [mm]'}        
        for k in self._keys:
             self._widgets[k] = ArrayInput(name=k, label=k, **inp[k])
        
        self._values['input'] = {}
        self._values['list'] = {}
        self._name_inp = ipy.Text(description='',value='', 
                                 layout=ipy.Layout(width='100px'))
        self._out = ipy.Output(layout=ipy.Layout(width='100%', 
                                                 border='1px solid')) 
    
    
    def _get_ori(self):
        """Retrieve geometry data from widgets."""
        vals = {}
        for key in self._widgets:
            vals[key] = self._widgets[key].value
        return vals
    
    def _on_add(self,b):
        txt = self._name_inp.value
#print(txt)
        if not txt:
            self.message('Define a unique name')
        elif not txt or txt in self._values['list']:
            self.message('The geometry name "{}" is already defined'.format(txt))
        else:
            self.add_geometry(txt, update_table=True)

    def _update_table(self):
        
        grid = self._create_list()
        self._out.clear_output()
        if len(self._values['list'].keys())>0:
            with self._out:
                display(grid)

    def _on_del_click(self, b):
        if b.value in self._values['list']:
            del self._values['list'][b.value]
            self._update_table()
            arg = {'list':self._values['list']}
            self._call_change(**arg)
    
    def _create_list(self):
        """Create a grid of widgets with info for all defined orientations."""
        grid_layout = ipy.Layout(grid_template_columns='50px 50px auto auto', width='30%')      
        hdr = ['','name','angles','scan']
        fmt = ['{}', '[{:g},{:g},{:g}]','[{:g},{:g},{:g}]']
        cell_layout = ipy.Layout(border='none')

        # header row
        hdr_inp = []
        for i in range(len(hdr)):
            L = ipy.HTML('<b>{}</b>'.format(hdr[i]),layout=cell_layout)
            hdr_inp.append(L)

        data = [hdr_inp]
        # collect rows of widgets for all defined orientations
        for key in self._values['list']:
            d = self._values['list'][key]
            rec = [key, d['angles'], d['scandir']]
            b = SButton(description='',icon='trash',value=key,
                            layout=ipy.Layout(width='30px'), tooltip='Remove')
            b.on_click(self._on_del_click)
            wrec = [b]
            for j in range(len(rec)):
                if isinstance(rec[j],np.ndarray):
                    v = list(rec[j])
                    L = ipy.HTML(fmt[j].format(*v,layout=cell_layout))
                else:
                    L = ipy.HTML(fmt[j].format(rec[j]),layout=cell_layout)
                wrec.append(L)
            data.append(wrec)

        flat_list = [item for sublist in data for item in sublist]
        gr = ipy.GridBox(flat_list, layout=grid_layout)        
        return gr
 
    def update_widgets(self, data):
        """Update input widgets from data.
        
        Parameters
        ----------
        data: dict
            Input data. The expected top-level keys are:
            
            input : str, optional
                Values for the geometry input widgets (four vectors).
            list : dict, optional
                Geometry datasets to be added in the table.
                Each dataset contains the four vectors defining a geometry.
        
        """
        if 'input' in data:
            dset = data['input']
            self._check_keys(dset)
            for key in self._keys:
                self._widgets[key].value = dset[key]
        
        if 'list' in data:
            lst = data['list']
            self._values['list'].clear()
            for key in lst:
                self._check_keys(lst[key])
                self.add_geometry(lst[key])
            self._update_table()
    
    def update_values(self):
        """Update input data from widgets.
        
        Does not change the list of defined orientations.

        """ 
        self._values['input'].clear()
        self._values['input'].update(self._get_ori())

    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        # header
        lbl = create_header('Sample orientation')
        lst = [lbl]
        # orientation data
        vals = list(self._widgets.values())
        for i in range(len(vals)):
            lst.append(vals[i].ui())
        # add item banner
        btn_add = ipy.Button(description='Add', 
                             layout=ipy.Layout(width='50px'))
        btn_add.on_click(self._on_add)
        hint = ipy.Label('Enter unique name (e.g. "radial", "axial", ...)')
        lst.append(ipy.HBox([hint, self._name_inp, btn_add]))
        # output area for table
        lst.append(self._out)
        super().show(widgets=lst, err=err, msg=msg)
        self._update_table()

    def add_geometry(self, name, update_table=False):
        """Add geometry using current input values."""
        if not name or name in self._values['list']:
            self.error('Provide a unique name for new orientation.')
            return
        vals = {}
        for key in self._widgets:
            vals[key] = self._widgets[key].value
        self._values['list'][name] =vals
        if update_table:
            self._update_table()
        arg = {'list':self._values['list']}
        self._call_change(**arg)
        self.error('')
        self.message('')
        
    def get_list(self):
        """Return the named list of added orientations.
        
        Attention
        ---------
        Returned is the reference to the internally stored dict, not 
        a copy. 
        
        Return
        ------
        dict
        """
        return self._values['list']


class UI_sampling(UI_base):
    """Input for a list of simulated sampling data sets.
    
    This UI permits to define a file with simulated sampling events and
    load required number of events from it. in addition, this UI 
    maintains a list of already loaded sampling files. 
    
    The values are returned by the method :meth:`get_values`. It returns
    named list with the keys:
        input:
            Filename, path and number of events for loading the sampling events . 
        list:
            A named list of sampling data sets added by the user.
            Only meta-data are included (file, path, nev). To get 
            corresponding list of Sampling objects, use the method
            :meth:`get_sampling` for a single sampling data set,
            or :meth:`get_list` for the named list of all 
            loaded data sets.
    
    Parameters
    ----------
    file : str
        Initial sampling file name
    path : str
        Initial directory name with sampling files
    nev : int
        Initial number of sampling points to be loaded.
    """
    
    style_lbl = {'description_width': '100px'}
    def __init__(self, file='', path='', nev=3000):
        super().__init__('sampling', ['file', 'path', 'nev'])        
        self._sampling = {}
        self._values['input'] = {}
        self._values['list'] = {}
        self._values['input']['file'] = file
        self._values['input']['path'] = path
        self._values['input']['nev'] = nev
        self._out = ipy.Output(layout=ipy.Layout(width='100%', 
                                                 border='1px solid')) 

        # file input
        fi = FileInput(name='file', 
                       file=self._values['input']['file'], 
                       path=self._values['input']['path'], 
                       label='File name',
                       tooltip='File with simulated sampling points.',
                       width_label=130)
        self._widgets['file'] = fi
        # num. of events
        #inp = create_input_int(name='file',label='Maximum events',
        #                       value=self._values['input']['nev'],
        #                       lmin=1000, lmax=100000, step=1000)   
        inp = ArrayInput(name='file',
                         label='Maximum events',
                         value=self._values['input']['nev'],
                         isInt=True, step=1000, lmin=3000,
                         width_label=130)
        self._widgets['nev'] = inp
        # ID string
        ninp = ipy.Text(description='Name',value='', 
                                 layout=ipy.Layout(width='200px'),
                                 tooltip = 'Provide unique name')
        self._widgets['name']  = ninp
    
    def _update_table(self):
        grid = self._create_list()
        self._out.clear_output()
        if len(self._values['list'].keys())>0:
            with self._out:
                display(grid)
        
    def _on_del_click(self, b):
        """Delete row from the list of loaded data sets.
        
        b.value should contain the key to the given record.
        """
        if b.value in self._values['list']:
            del self._values['list'][b.value]
        if b.value in self._sampling:
            del self._sampling[b.value]
        self._update_table()
        self.error('')
        arg = {'list':self._values['list']}
        self._call_change(**arg)
    
    def _on_load_click(self, b):
        """Select file and load sampling."""
        self.update_values()
        name = self._widgets['name'].value
        self.add_sampling(name, self._values['input'], update_table=True)
        
    def _create_list(self):
        """Create a grid of widgets with info for loaded sampling datasets."""
        grid_fmt = '50px 80px auto repeat(7,70px)'
        grid_layout = ipy.Layout(grid_template_columns=grid_fmt)      
        hdr = ['','ID','file','nrec','&lambda;','2theta','dhkl', 
               'width_x','width_y','width_z']
        fmt = ['{}', '{}', '{:d}'] + 3*['{:.5g}'] + 3*['{:.2f}']
        cell_layout = ipy.Layout(border='none')

        # header row
        hdr_inp = []
        for i in range(len(hdr)):
            L = ipy.HTML('<b>{}</b>'.format(hdr[i]),layout=cell_layout)
            hdr_inp.append(L)

        data = [hdr_inp]
        # collect rows of widgets for loaded sampling data sets
        for key in self._values['list']:
            #d = self._values['list'][key]
            s = self._sampling[key]
            w = s.src['width']
            rec = [key, s.file, s.src['nrec'], s.src['wav'], s.src['tth'],
                   s.src['dmean'], w[0], w[1], w[2]]
            b = SButton(description='', icon='trash', value=key,
                            layout=ipy.Layout(width='30px'), 
                            tooltip='Remove')
            b.on_click(self._on_del_click)
            wrec = [b]
            for j in range(len(rec)):
                L = ipy.HTML(fmt[j].format(rec[j]),layout=cell_layout)
                wrec.append(L)
            data.append(wrec)

        flat_list = [item for sublist in data for item in sublist]
        gr = ipy.GridBox(flat_list, layout=grid_layout)        
        return gr
       
    def update_widgets(self, data):
        """Update widgets from data.
        
        If data has the key 'input', the input widgets are updated.
        
        If data hse the key 'list', the list of sampling data sets is updated
        and corresponding data loaded.
        
        `input` contains keys [`file`,`path`,`nev`] for file and path names
        and number of events to be loaded.
        
        `list` is a named list (dict) with inputs for sampling data
        to be loaded.

        Parameters
        ----------
        data : dict
            Expected keys are [`input`,`list`] (both optional)
        """
        if 'input' in data:
            self._widgets['file'].value = data['input']
            self._widgets['nev'].value = data['input']['nev']
     
        if 'list' in data:
            lst = data['list']
            self._values['list'].clear()
            self._sampling.clear()
            for key in lst:
                self.add_sampling(key, lst[key], update_table=False)
            self._update_table()    
    
    def update_values(self):
        """Update input data from widgets.
        
        Does not change the list of defined sampling sets.

        """ 
        file = self._widgets['file'].value
        self._values['input']['file'] = file['file']
        self._values['input']['path'] = file['path']
        self._values['input']['nev'] = self._widgets['nev'].value

    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        lst = []
        lbl = create_header('Sampling points')      
        # load button
        btn = SButton(description='Load',  icon='file', 
                          layout=ipy.Layout(width='10%'))     
        btn.on_click(self._on_load_click)
        
        # collect boxes and display
        box_add = ipy.HBox([self._widgets['nev'].ui(), 
                            self._widgets['name'], 
                            btn ])
        box = ipy.VBox([lbl, self._widgets['file'].ui(), box_add])
        lst.append(box)
        lst.append(create_header('Loaded sampling data',size='-1'))  
        lst.append(self._out)
        super().show(widgets=lst, err=err, msg=msg)
        # display updated table 
        self._update_table()
    
    def add_as(self, name, update_table=False):
        """Load and add sampling on the list according to the current input."""
        self.add_sampling(name, self.get_values()['input'], 
                          update_table=update_table )
    
    def add_sampling(self, name, data, update_table=False):
        """Load sampling and add it to the list.
        
        Parameters
        ----------
        name : str
            Unique sampling name
        data : dict
            file, path and nev parameters. `data` is passed as keyward aguments
            to :meth:`stressfit.commands.load_sampling`.
        update_table: bool
            if true, replot the list of loaded sampling sets
        """   
        if not name or name in self._values['list']:
            self.error('Provide a unique name for new sampling.')
            return
        
        self.error('')
        try:
            self._check_keys(data)
            sampling = comm.load_sampling(**data)
            if sampling is not None:
                data['file'] = sampling.src['file']
                data['path'] = sampling.src['path']
                data['nev'] = sampling.src['nrec']
                #inp = {k: data[k] for k in data.keys() & {'file', 'path'}}
                self._widgets['file'].value = data
                self._values['input'].update(data)
                self._values['list'][name]=data.copy()
                self._sampling[name]=sampling
                if update_table:
                    self._update_table()  
                arg = {'list':self._values['list']}
                self._call_change(**arg)
        except Exception as e:
            self.error(e)
            
    def get_sampling(self, key):
        """
        Get selected sampling data set.

        Parameters
        ----------
        key : str
            Key for the seleted data set.

        Returns
        -------
        Sampling
            An instance of the Sampling class.

        """
        if key in self._sampling:
            return self._sampling[key]
        else:
            return None
        
    def get_list(self):
        """Return the named list of loaded sampling data sets.
        
        Attention
        ---------
        Returned is the reference to the internally stored dict, not 
        a copy. 
        
        Return
        ------
        dict
        """
        return self._sampling


class UI_attenuation(UI_base):
    """Dialog for beam attenuation.

    Parameters
    ----------
    workspace: stressfit.dataio.Workspace
        Reference to the workspace object.
    data: dict
        Initial input values (optional)
    
    """
    
    _att_types = ['value', 'table']
    def __init__(self, workspace, data=None):
        super().__init__('attenuation', ['type','value','table']) 
        self.wk = workspace
        self._values['table'] = 'Fe_mu.dat'
        self._values['value'] = 1.1
        self._values['type'] = UI_attenuation._att_types[1]

        radio = SRadioButtons(name='type', 
                             options=UI_attenuation._att_types,
                             value=self._values['type'],
                             description='')
                
        #val = create_input_float(name='value', 
        #                         label='att. coefficient [1/cm]',
        #                         value=self._values['value'], 
        #                         lmin=0.0)        
        
        val = ArrayInput(name='value', 
                   label='Value',
                   value=self._values['value'], 
                   hint='attenuation coefficient [1/cm]',
                   width_label=80)
        
        fi = FileInput(name='table',
                       file=self._values['table'], 
                       path=self.wk.full_path('tables', as_posix=True), 
                       label='Table',
                       tooltip='Lookup table for attenuation coefficient [1/cm] vs. wavelength [AA].',
                       width_label=80)
        
        radio.observe(self._type_change)
        self._widgets[radio.name] = radio
        self._widgets[val.name] = val
        self._widgets[fi.name] = fi
        if data is not None:
            self.update_widgets(data)
            self.update_values()

    def _type_change(self, change):
        if (change['name']=='value'):
            t = change['new']
            is_value = (t=='value')
            self._widgets['table'].disabled = is_value
            self._widgets['value'].disabled = not is_value
     
            
    def update_widgets(self, data):
        """Update widgets from data.
        
        Expected input keys are:
        
        value : float
            Attenuation coefficient [1/cm]
        table : str
            File with a table of att. coefficients:
                (wavelength [AA], value [1/cm])
        type : str
            'value' to use the scalar value, 'table' to use the lookup table.

        Parameters
        ----------
        data : dict
            Input data.
        """
        for key in data:
            if key in self._widgets:
                self._widgets[key].value = data[key]
        
    def update_values(self):
        """Update input data from widgets.
        
        Does not change the list of defined sampling sets.

        """ 
        for key in self._widgets:
            self._values[key] = self._widgets[key].value
        
        
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        lbl1 = create_header('Attenuation')
        box = [lbl1] 
        box.append(self._widgets['type'])
        box.append(self._widgets['value'].ui())
        box.append(self._widgets['table'].ui(width='80%'))
        super().show(widgets=box, err=err, msg=msg)
        self._type_change({'name':'value', 'new':self._values['type']})
        
    def get_attenuation(self):
        """Return attenuation input.
        
        Returns
        -------
        float or str
            Either single value [1/cm] or the file name.
        """
        self.update_values()
        
        if self._values['type'] == 'value':
            att = self._values['value']
        else:
            val = self._widgets['table'].value
            file = val['file']
            path = val['path']
            att = dataio.load_data(file, kind='tables', path=path)
        return att

class UI_plot_scene(UI_base):
    """Plot scene with sample geometry and sampling events.
    
    A dialog for plotting of the sample contours with scattering geometry 
    and sampling events in several projections. 
    
    Parameters
    ----------
    samplings: dict
        List of loaded sampling data sets.
        as defined by the class :class:`UI_sampling` (use :meth:`get_values`)
    geometries: dict
        List of oriantations, as defined by the class :class:`UI_orientations`
        (use :meth:`get_values`)
    rang: int
        Display range in mm. 
    proj: int
        Initially selected projection (0 ..2) 
    nev: int
        Number of sampling events to show.
        
    Example
    -------
    Provided that the objects `s` and `g` are the dialogs for sampling
    and geometry, respectively, than create and display this dialog
    as 
    
    `plot = UI_plot_scene(s.get_list(), g.get_list())`
    
    `plot.show()`
    
    
    """
    
    def __init__(self, samplings, geometries, rang=16, proj=1, nev=3000):
        super().__init__('scene', ['sampling','ori','nrec','proj','range'])    
        self.header = create_header('Plot scene')
        # sampling selection
        wdg = create_select(name='sampling', label='Sampling', 
                            options=[], 
                            width_label=100, width_drop=100)
        self._widgets['sampling'] = wdg
        # oreiantation selection
        wdg =  create_select(name='ori', label='Orientation', 
                             options=[], 
                             width_label=100, width_drop=100)
        self._widgets['ori'] = wdg
        # number of events to plot
        wdg = create_input_int(name='nrec', label='Events',
                               value=nev, lmin=1000, lmax=100000, step=1000,
                               width_label=100, width_num=100)
        self._widgets['nrec'] = wdg 
        # projection plane
        wdg = create_select(label='Projection: ', 
                            options=[('z,y',0),('x,z',1),('x,y',2)],
                            value = proj,
                            width_label=100, width_drop=100)
        self._widgets['proj'] = wdg 
        # plot range
        wdg = ipy.IntSlider(value=rang, min=3, max=50, step=1, 
                            description='Range: ')
        self._widgets['range'] = wdg 
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none')) 
        
        # select options
        self.options_sampling = []
        self.options_ori = []
        self.update_sampling_options(samplings)
        self.update_geometry_options(geometries)
        
    def update_sampling_options(self, samplings):
        """Update the selection list of sampling data."""
        self.samplings = samplings
        self.options_sampling.clear()
        for key in self.samplings:
            self.options_sampling.append((key,key))
        self._widgets['sampling'].options = self.options_sampling
        self._widgets['sampling'].index = 0
            
    def update_geometry_options(self, geometries):
        """Update the selection list of geometries."""
        self.geometries = geometries
        self.options_ori.clear()
        for key in self.geometries:
            self.options_ori.append((key,key))
            self._widgets['ori'].options = self.options_ori
            self._widgets['ori'].index = 0

    def _on_replot(self,b):
        # Plot experiment geometry 
        # (red arrow shows motion of sampling points in stationary sample)
        
        # read selected geometry
        geometry = self.geometries[self._widgets['ori'].value]
        sampling = self.samplings[self._widgets['sampling'].value]
        self._call_init(sampling=sampling, geometry=geometry)


       #  set_input(geometry, sampling)
        self._out.clear_output()
        with self._out:
            comm.plot_scene(self._widgets['nrec'].value, 
                            filename='', 
                            rang=2*[self._widgets['range'].value],
                            proj=self._widgets['proj'].value)    
        
    def update_widgets(self, data):
        """Update widgets from input data.""" 
        self._update_options()
        if 'sampling' in data:
            self.update_sampling_options(self, data['sampling']) 
        elif 'geometry' in data:
            self.update_geometry_options(self, data['geometry'])  
        else:
            for key in self._widgets:
                if key in data:
                    self._widgets[key].value = data[key]

    def update_values(self):
        """Update input data from widgets.""" 
        for key in self._widgets:
            self._values[key] = self._widgets[key].value
            
    def show(self, err=None, msg=None, **kwargs):
        """Display the input collection."""
        lbl = create_header('Plot scene')
        btn_replot = ipy.Button(description='Replot',
                                layout=ipy.Layout(width='80px'))
        btn_replot.on_click(self._on_replot)
        box1 = ipy.VBox([self._widgets['sampling'], 
                         self._widgets['ori'], 
                         self._widgets['proj'],
                         self._widgets['nrec']])
        box_plot = ipy.HBox([self._widgets['range'], btn_replot])
        box2 = ipy.VBox([box_plot, self._out])
        lst = [lbl, ipy.HBox([box1, box2])]
        super().show(widgets=lst, err=err, msg=msg)
 
 
#%% Top level UI class

class UI():
    """Top level class generating the Jupyter notebook interface."""
    
    _registered_ui = ['workspace',
                      'shape',
                      'geometry',
                      'sampling',
                      'attenuation',
                      'scene']
    def __init__(self):
        self.setup = {}
        self.ui = {}
        self.wk = dataio.workspace()
        self._msg = ipy.Output(layout=ipy.Layout(width='100%'))
        self._err = ipy.Output(layout=ipy.Layout(width='100%'))
        # setup dialog
        self._add_input_ui(UI_workspace())
        self._add_input_ui(UI_shape(select=shapes.Tube))
        self._add_input_ui(UI_geometry())
        self._add_input_ui(UI_sampling())
        self._add_input_ui(UI_attenuation(self.wk))
        
        # initialize
        self.ui['geometry'].add_geometry('radial', update_table=False)
        inpdir = self.ui['workspace'].get_values()['data']
        sampl = {'path':inpdir, 'file':'events_S_1mm.dat', 'nev':3000}
        self.ui['sampling'].add_sampling('1x1x5mm',sampl, update_table=False )
        
        # execution dialogs
        self.exec = {}
        self.ui['scene'] = UI_plot_scene(
            self.ui['sampling'].get_list(),
            self.ui['geometry'].get_list())
        
        # set event handlers
        self.ui['workspace'].set_on_change(self._change)
        self.ui['shape'].set_on_change(self._change)
        self.ui['geometry'].set_on_change(self._change)
        self.ui['sampling'].set_on_change(self._change)
        
        self.ui['scene'].set_on_init(self._init)
    
    def _add_input_ui(self, ui):
        if ui.name in UI._registered_ui:
            self.ui[ui.name] = ui
            self.setup[ui.name] = ui.get_values()
        else:
            msg = 'Attempt to add non-registered ui: {}'
            raise Exception(msg.format(ui.name))
        
    def display(self):
        """Display all input blocks."""
        # output tabs layout
        layout=ipy.Layout(width='100%')        
        tabs_data = {}
        tabs_data['workspace'] = {'title':'Workspace',
                                  'ui':[self.ui['workspace']]}
        tabs_data['geometry'] = {'title':'Geometry',
                                  'ui':[self.ui['shape'],
                                        self.ui['geometry']]}
        tabs_data['sampling'] = {'title':'Sampling',
                                  'ui':[self.ui['sampling'],
                                        self.ui['scene']]}
        tabs_data['material'] = {'title':'Material',
                                  'ui':[self.ui['attenuation']]}
        # create tab container
        tab = ipy.Tab() 
        
        # create output tabs
        tabs = {}
        for key in tabs_data:
            tabs[key] = ipy.Output(layout=layout)
        tab.children = list(tabs.values())
        
        # set tab titles
        keys = list(tabs_data.keys())
        for key in keys:
            i = keys.index(key)
            tab.set_title(i, tabs_data[key]['title'])
            
                
        # display all
        display(tab)
        keys = list(tabs_data.keys())
        for key in keys:
            for ui in tabs_data[key]['ui']:
                with tabs[key]:
                    ui.show(msg=self._msg, err=self._err)
        display(self._err)
        display(self._msg) 
    
    def _change(self, obj, **kwargs):
        data = self.setup[obj.name]
        if obj.name == 'shape':
            if 'shape' in kwargs:
                data.update(**obj.get_values())
                comp = self.ui['shape'].create_shape()
                comm.set_shape(comp)
            else:
                data['param'].update(**kwargs)
                comm.set_shape(None,**kwargs)
        elif obj.name == 'workspace':
            if 'work' in kwargs:
                data.update(**obj.get_values())
            else:
                data.update(**kwargs)
        elif obj.name in ['geometry']:
                data.update(**kwargs)
                new_options = self.ui['geometry'].get_list()
                self.ui['scene'].update_geometry_options(new_options)
        elif obj.name in ['sampling']:
                data.update(**kwargs)
                new_options = self.ui['sampling'].get_list()
                self.ui['scene'].update_sampling_options(new_options)

    def _init(self, obj, **kwargs):
        """Initialize stressfit.
        
        Callback used by IU groups to initialize stressfit.
        """        
        if obj.name == 'scene':
            if 'sampling' in kwargs:
                comm.set_sampling(kwargs['sampling'])
            if 'geometry' in kwargs:
                comm.set_geometry(kwargs['geometry'])
            

    def save(self, filename=''):
        """Save input data in JSON format."""
        out = {'ui': self.setup} 
        txt = json.dumps(out,indent=4)
        if filename:
            wkp = self.wk.get_paths(keys=['work'])
            file = wkp['work'].joinpath(filename)
            f = open(file, 'w')
            f.write(txt)
            f.close()
        else:
            with self._msg:
                print(txt)

        
