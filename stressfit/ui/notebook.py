"""Classes and functions for creating notebook interface.

Created on Tue Oct  5 11:38:15 2021
@author: Jan Saroun, saroun@ujf.cas.cz
"""

# TODO implement _create_widgets for all descendants of UI_base

# TODO exp. data input
# TODO model definition - intensity
# TODO model definition - strain
# TODO fit intensity
# TODO fit strain
# TODO input - compliance
# TODO model definition - stress
# TODO fit stress
# 

#import abc
import abc
import ipywidgets as ipy
import json
import copy
from IPython.display import display
from .widgets import DirInput, FileInput, ArrayInput, SelectInput
from .widgets import create_header, create_select, create_input_int 
from .widgets import create_input_float, create_checkbox, create_text
from .widgets import SButton, SRadioButtons, SCheckbox
from .widgets import choose_file, choose_file_save
from pathlib import Path as _Path
import numpy as np
import stressfit.commands as comm
import stressfit.shapes as shapes
import stressfit.dataio as dataio
import stressfit.graphs as gr



def _has_keys(args, keys):
    """Verify that all args are in keys."""
    return all (k in args for k in keys)

#%% Base abstract classes for input collections


class UI_base:
    """Common parent class for all input collection.
    
    Parameters
    ----------
    name : str
        Unique name for the instance.
    keys : list
        List of keys for allowed values and widgets.
    kwargs : dict
        Arguments passed to :meth:`_init_values`.
    
    Abstract methods
    ----------------
    _init_values
        Initiate values, called by constructor.
    _create_widgets
        Create value widgets, called by :meth:`show`.
    _create_ui
        Create a list of all ipywidgets to be placed in the main VBox.

    """
    
    def _err_ni(src):
        msg = "{}(): Subclass must implement abstract method".format(src)
        raise NotImplementedError(msg)
    
    def __init__(self, name, keys, **kwargs):   
        # define layout parameters - subclasses should override it
        self._err = None # error output, can be defined by ui()
        self._msg = None # message output, can be defined by ui()
        self._debug = False
        self.uiparam = {}
        self.uiparam['title'] = 'Title'
        self.uiparam['title_size'] = '+1'
        self.uiparam['title_color'] = '#35446B'
        self._name = name # unique instance name
        self._keys = keys # list of allowed value keys
        self._widgets = {} # dict of input widgets
        self._values = {} # dict of widget values
        self._buttons = [] # Command buttons next to the collection label 
        self._on_change = None
        self._on_init = None
        self._options = {} # container for selection parmeters
        self._init_values(**kwargs)
        
# abstract methods

    @abc.abstractmethod
    def _init_values(self, **kwargs):
        """Initiate values.
        
        Called by constructor. 
        """
        UI_base._err_ni('_init_values')

    @abc.abstractmethod
    def _create_widgets(self, **kwargs):
        """Create widgets registered in self._widgets dict for all input keys.
        
        Input keys are defined by the constructor argument 'keys'.
        """
        UI_base._err_ni('_create_widgets')

    @abc.abstractmethod
    def _create_ui(self, **kwargs):
        """Create a top-level ipywidget to be placed in the main VBox.
        
        Envelop the input widgets created by :meth:`_create_widgets` 
        and add other widgets like titles control buttons etc. 
                
        Returns
        -------
        widget
            A top level widget enveloping all inputs, title, buttons, etc.
        
        """
        UI_base._err_ni('_create_ui')

# private methods

    def _dbg_message(self, txt, clear=False):
        if self._debug and self._msg:
            if clear:
                self._msg.clear_output()
            with self._msg:
                print(txt)

    def _check_keys(self,data):
        """Verify that data contains all required keys."""
        if not _has_keys(data,self._keys):
            msg = 'Required key(s) missing in the input data: {}'
            raise Exception(msg.format(self._keys))
 
    def _call_change(self, **kwargs):
        """Call when the UI changes state."""
        if self._on_change:
            self._on_change(self, **kwargs)
    
    def _call_init(self, **kwargs):
        """Call on the UI initialization."""
        if self._on_init:
            self._on_init(self, **kwargs)

# properties
        
    @property
    def name(self):
        """ID string for the input block."""
        return self._name


# public methods

    def notify(self):
        """Raise notification after change.
        
        UI_base does nothing. Subclasses may need to implement its own actions.
        """
        

    def set_values(self, values):
        """Set values.
        
        Set internally stored data in `self._values`. 
        To update widgets for the new values, call :meth:`update_widgets`.
        
        Parameters
        ----------
        values : dict
            New values.
        """
        for key in self._keys:
            if key in values:
                self._values[key] = copy.deepcopy(values[key])

    def get_values(self):
        """Return actual input values as dict."""
        self.update_values()
        return self._values
    
    def update_widgets(self):
        """Update input widgets from values."""        
        for key in self._widgets:
            if key in self._values:
                self._widgets[key].value = self._values[key]

    def update_values(self):
        """Update values from widgets."""        
        for key in self._widgets:
            self._values[key] = copy.deepcopy(self._widgets[key].value)

    def update_options(self, name, options):
        """Update the selection list."""
        self._dbg_message('update_options: name={}'.format(name))
        self._dbg_message('{}'.format(options))
        if not name in self._options:
            return
        if not name in self._widgets:
            msg = 'No selection widget corresponding to: {}'
            raise Exception(msg.format(name))    
        # avoid clear if update options is called on itself 
        if options != self._options[name]:
            self._options[name].clear()
            self._options[name].update(options)
            self._dbg_message('updated')
        
        # construct selection widget options
        select_options = []
        for key in options:
            select_options.append((key,key))
            
        # update widget    
        val = self._widgets[name].value   
        self._widgets[name].options = select_options
        self._dbg_message('{}'.format(self._widgets[name].options))
        # try to restore previous selection
        
        try:
            if val in self._widgets[name].options:
                self._widgets[name].value = val
            elif len(self._widgets[name].options)>0:
                self._widgets[name].index = 0
        except Exception as e:
            print(e)

    def message(self, txt, clear=True):
        """Print message in the message output, if defined."""
        if self._msg:
            if clear:
                self._msg.clear_output()
            with self._msg:
                print(txt)
        else:
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
        else:
            print(txt)

    def set_on_change(self, on_change):
        """Set callback to be executed when the UI changes state."""
        self._on_change = on_change
        
    def set_on_init(self, _on_init):
        """Set callback to be executed when the UI needs to initialize."""
        self._on_init = _on_init
              
    def show(self, err=None, msg=None, **kwargs):
        """Display VBox container with all input widgets."""        
        self._err = err
        self._msg = msg
        self._create_widgets(**kwargs)
        ui = self._create_ui(**kwargs)
        hdr = create_header(self.uiparam['title'],
                      color=self.uiparam['title_color'],
                      size=self.uiparam['title_size'])
        top = [hdr] + self._buttons
        top_wdg = ipy.HBox(top)
        layout = ipy.Layout(margin='0px 0px 20px 0px')
        display(ipy.VBox([top_wdg, ui], layout=layout))


class UI_base_list(UI_base):
    """Extends UI_base by adding a list of input data sets."""
    
    def __init__(self, name, keys, 
                 list_template='',
                 list_hdr='',
                 list_fmt='',
                 **kwargs):
        super().__init__(name, keys, **kwargs)
        # define additional layout parameters       
        self.uiparam['list_title'] = 'List title'
        self.uiparam['list_title_size'] = '-1'
        self.uiparam['list_border'] = '1px solid'
        self.uiparam['list_width'] = '100%'
        self.uiparam['add_button_label'] = 'Add'
        self.uiparam['add_button_icon'] = None
        self.uiparam['add_button_width'] = '50px'
        self.uiparam['add_name_width'] = '200px'
        self.uiparam['add_label'] = 'Enter unique name'
        self.uiparam['add_label_width'] = '200px'
        
        # save format info for the table of data sets.
        self._template = list_template
        self._hdr = list_hdr
        self._fmt = list_fmt



    def _init_values(self, **kwargs):
        """Initiate values within constructor.
         
        Only creates dict `_values['list']`, `_values['input']` and `_data`.
        Subclasses should implement actual data initialization.
        """
        #print('UI_base_list._init_values\n{}'.format(kwargs))
        self._values['list'] = {}
        self._values['input'] = {}
        self._data = {}
        
        
    def _is_unique(self, name):
        """Check that name is a unique key for `self._values['list']`."""
        out = False
        if not name:
            self.error('Define a unique name')
        elif name in self._values['list']:
            self.error('The item ID "{}" is already defined'.format(name))
        else:
            out = True
        return out
    
    def _update_table(self):
        """Redraw the table of data sets in the `self._out` area."""
        grid = self._create_list()
        try:       
            self._out.clear_output()
            if len(self._values['list'].keys())>0:
                    with self._out:
                        display(grid)
        except:
            pass

    def _delete(self, key):
        """Delete given item from the list.
        
        If successful, updates the table, but does not call `on_change`.
        """
        if key in self._values['list']:
            del self._values['list'][key]
            if key in self._data:
                del self._data[key]
            self._update_table()
            return True
        else:
            return False
            
    def notify(self):
        """Raise notification after change.
        
        UI_base_list calls `self._call_change(list=self._values['list'])`.
        """
        arg = {'list':self._values['list']}
        self._call_change(**arg)
        
    def _on_del_click(self, b):
        """Delete the table row and call `on_change`."""
        if self._delete(b.value):
            self.notify()

    def _on_add(self,b):
        """Call :meth:`add_data`.

        `on_change` is sent at the end with 'list' and 'data' arguments.
        """
        txt = self._name_inp.value
        
        self._dbg_message('Add {}'.format(txt), clear=True)
        
        self.add_data(txt, update_table=True)
        self.notify()

    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given ID key.
        
        By default, return a list with key followed by all items from 
        `_values['list'][key]`. Subclasses may need to modify this behaviour.
        
        """
        d = self._values['list'][key]
        rec = [key]
        for k in d.keys():
            rec.append(d[k])
        return rec        

    def _create_list(self):
        """Create a grid of widgets with info for the input list."""
        # grid for delete_btn and all items returned by _get_display_record
        grid_layout = ipy.Layout(grid_template_columns='50px ' + self._template)      
        hdr = [''] + self._hdr
        fmt = self._fmt
        #print(fmt)
        cell_layout = ipy.Layout(border='none')

        # header row
        hdr_inp = []
        for i in range(len(hdr)):
            L = ipy.HTML('<b>{}</b>'.format(hdr[i]),layout=cell_layout)
            hdr_inp.append(L)

        data = [hdr_inp]
        # collect rows of widgets for all defined orientations
        self._dbg_message('_create_list')
        for key in self._values['list']:
            rec = self._get_display_record(key)
            
            self._dbg_message('\t{}'.format(rec))
            
            b = SButton(description='',icon='trash',value=key,
                            layout=ipy.Layout(width='30px'), tooltip='Remove')
            b.on_click(self._on_del_click)
            wrec = [b]
            for j in range(len(rec)):
                if isinstance(rec[j],list):
                    L = ipy.HTML(fmt[j].format(*rec[j],layout=cell_layout))
                else:
                    L = ipy.HTML(fmt[j].format(rec[j]),layout=cell_layout)
                wrec.append(L)
            data.append(wrec)

        flat_list = [item for sublist in data for item in sublist]
        gr = ipy.GridBox(flat_list, layout=grid_layout)        
        return gr

    def set_values(self, values):
        """Set values.
        
        Set internally stored data in `self._values` and `self._data`. 
        To update widgets for the new values, call :meth:`update_widgets`.
        
        Parameters
        ----------
        values : dict
            Input data. The expected top-level keys are:
            
            input : str, optional
                Keys and values for the input widgets.
            list : dict, optional
                Data sets to be added in the table. 
        """
        if 'input' in values:
            for key in self._keys:
                if key in values['input']:
                    self._values['input'][key] = copy.deepcopy(values['input'][key])
                    #print('set_values {}\n\t{}\n\t{}'.format(key,self._values['input'][key],values['input'][key]))
        if 'list' in values: 
            lst = values['list']
            self._values['list'].clear()
            self._data.clear()
            for key in lst:
                self.add_data(key, data=lst[key])                
            
    def update_values(self):
        """Update input data from widgets."""
        # self._values['input'].clear()
        for key in self._keys:
            if key in self._widgets:
                self._values['input'][key] = copy.deepcopy(self._widgets[key].value)
    
    def update_widgets(self):
        """Update input widgets from data."""
        for key in self._keys:
            if key in self._widgets:
                self._widgets[key].value = self._values['input'][key]
        self._update_table()

    def add_data(self, name, data=None, update_table=False):
        """Add a dataset to the list.
        
        By default, put all widgets values on the list.
        Subclasses may need to override this e.g. by loading data files etc.
        
        Parameters
        ----------
        name : str
            ID string for new item. It must be a unique ID string.
        data : dict
            New data set to be added to the table.
            If None, then values are taken from corresponding widgets.
        update_table : bool
            If True, call :meth:`_update_table` to redraw the table content.
        """
        if not self._is_unique(name):
            return
        if data is None:
            self.update_values()
        else:
            try:
                self._check_keys(data)
            except Exception as e:
                print(e)
            self.set_values({'input':data})
        self._values['list'][name] = copy.deepcopy(self._values['input'])
        
        self._dbg_message('UI_base_list.add_data: {}'.format(name))
        self._dbg_message('{}'.format(self._values['list'][name]))
        if update_table:
            self._update_table()
        
        
    def get_data(self, key=None):
        """
        Get selected data set.

        Parameters
        ----------
        key : str
            Key for the seleted data set.

        Returns
        -------
        Object
            Loaded data set for given key.
            if key==None then return all ._data

        """
        if key is None:
            return self._data
        elif key in self._data:
            return self._data[key]
        else:
            return None
        
    def get_list(self):
        """Return the named list of loaded data sets.
        
        Attention
        ---------
        Returned is the reference to the internally stored dict, not 
        a copy. 
        
        Returns
        -------
        dict
        """
        return self._values['list']

    def show(self, err=None, msg=None, **kwargs):
        """Display VBox container with all input widgets.
        
        In addition to UI_base, display the table with the list
        of data sets, text field and a button for adding new items.
        
        Parameters
        ----------
        err : ipywidgets.Output
            Output area for errors provided by UI application
        msg : ipywidgets.Output
            Output area for messages provided by UI application
        
        kwargs : dict
            Parameters passed to :meth:`_create_ui`
        """               
        self._err = err
        self._msg = msg
        
        # create value widgets and place them in vertical list
        self._create_widgets(**kwargs)
        ui = self._create_ui(**kwargs)
        
        # create title
        hdr = create_header(self.uiparam['title'],
                      color=self.uiparam['title_color'],
                      size=self.uiparam['title_size'])

        # create add button
        btn_lbl = self.uiparam['add_button_label']
        btn_w = self.uiparam['add_button_width']
        btn_add = ipy.Button(description=btn_lbl, 
                             layout=ipy.Layout(width=btn_w))
        btn_add.on_click(self._on_add)
        # create add label
        hint_w = self.uiparam['add_label_width']
        hint = ipy.Label(self.uiparam['add_label'],
                         layout=ipy.Layout(width=hint_w))
        # create name input
        name_w = self.uiparam['add_name_width']
        self._name_inp = ipy.Text(description='',value='', 
                                 layout=ipy.Layout(width=name_w))
        
        # create list header
        list_hdr = create_header(self.uiparam['list_title'],
                                 size=self.uiparam['list_title_size'],
                                 color=self.uiparam['title_color'])
        
        # create table output area
        tab_border = self.uiparam['list_border']
        tab_width = self.uiparam['list_width']
        self._out = ipy.Output(layout=ipy.Layout(width=tab_width, 
                                                 border=tab_border))
        
        # make HBox with add button and name field
        add_box = ipy.HBox([hint, self._name_inp, btn_add])
        
        # display everything
        top = [hdr] + self._buttons
        top_wdg = ipy.HBox(top)
        layout = ipy.Layout(margin='0px 0px 20px 0px')
        display(ipy.VBox([top_wdg, ui, add_box, list_hdr, self._out], 
                         layout=layout))
        # update table
        self._update_table()      
        
#%% Input collections
        
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
    
    def __init__(self, name, exclude = ['instruments']):
        # NOTE: the keys list must match dataio.__path_keys
        self.wk = dataio.workspace()
        super().__init__(name, dataio.Workspace.types, exclude=exclude)
        self.uiparam['title'] = 'Workspace'
        self._out = ipy.Output(layout=ipy.Layout(width='100%'))     
        
    def _init_values(self, **kwargs):
        self._excluded_keys = kwargs['exclude']
        self._values.update(self.wk.get_paths())
    
    def _create_widgets(self, **kwargs):
        
        for key in self._keys:
            if not key in self._excluded_keys:
                [lbl, ttp] = UI_workspace._inp[key]
                inp = DirInput(name=key, 
                               path=self._values[key], 
                               label=lbl, 
                               tooltip=ttp)
                inp.observe(self._on_path_change)
                self._widgets[key] = inp
        self._widgets['work'].txt.disabled=True
        self._widgets['work'].callback = self._on_workspace_change       
    
    def _create_ui(self, **kwargs):    
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
        
        layout = ipy.Layout(display='flex',flex_flow='row',
                                border='none',width='100%')
        appl_box = ipy.HBox([appl, save, load, rst], layout=layout)
        
        lst = []
        for w in self._widgets:
            lst.append(self._widgets[w].ui())
        lst.append(appl_box)
        lst.append(self._out)
        return ipy.VBox(lst)
    
    def notify(self):
        """Raise notification after change.
        
        Calls `self._call_change(work=self.wk.get_paths())`.
        """
        data = self.wk.get_paths()
        self._call_change(**{'work':data['work']})
    
    def _on_workspace_change(self,s):
        with self._out:
            res = self.wk.change_workspace(self._values['work'])
            data = self.wk.get_paths()
            self.update_widgets(data)
            self.update_values()
            if res:
                self.message('Workspace configuration reloaded.')
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
                self.message(msg.format(wkp['work'],dataio.Workspace.cfg_name))
        except Exception as e:
            print(e)
    
    def _on_load(self, b):
        """Re-load workspace configuration from workspace root directory."""
        with self._out:
            if self.wk.load():
                data = self.wk.get_paths()
                self.update_widgets(data)
                self.update_values()
                self.message('Workspace configuration reloaded.')
            else:
                wkp = self.wk.get_paths(keys=['work'])
                msg = 'Workspace configuration file {} not found in {}'
                self.error(msg.format(dataio.Workspace.cfg_name, wkp['work']))     
       
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


class UI_options(UI_base):
    """UI for setting program options."""    
    
    def __init__(self, name, prefix='', save=True):
        super().__init__(name, ['prefix', 'save'],
                         prefix=prefix, save=save)         
        self.uiparam['title'] = 'Options'
         
    def _init_values(self, **kwargs):
        self._values['prefix'] = kwargs['prefix']
        self._values['save'] = kwargs['save']
        
    def _create_widgets(self, **kwargs):        
        # auto-save
        wdg = create_checkbox(name='save', label='Save',
                              value=self._values['save'],
                              width_label=100, width_box=20,
                              tooltip='Auto-save results',
                              indent=False)
        self._widgets['save'] = wdg
        wdg.observe(self._observe,'value', type='change')
        # prefix for output file names
        wdg = create_text(name='prefix', label='Prefix', 
                          value=self._values['prefix'],
                          width_label=60, width_box=100,
                          tooltip='Prefix for output file names')
        self._widgets['prefix'] = wdg
        wdg.observe(self._observe,'value', type='change')
    
    def _create_ui(self, **kwargs):
        line1 = ipy.HBox([self._widgets['save'], self._widgets['prefix']])
        return line1
    
    def _observe(self, change):
        if change['name']=='value':
            owner = change['owner']
            self._values[owner.name] = change['new']


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
    
    def __init__(self, name, select=shapes.Plate):
        self.select = select
        self.shapes = dataio.load_config('shapes')
        super().__init__(name, ['shape', 'param'])        
        self.uiparam['title'] = 'Sample shape'
        self._out = ipy.Output(layout=ipy.Layout(border='1px solid'))

    def _init_values(self, **kwargs):
        self._values['shape'] = self.select
        self._values['param'] = self.shapes[self.select]['param']
               
    def _create_widgets(self, **kwargs):
        options = []       
        for key in self.shapes:
            options.append((self.shapes[key]['name'],key))
        
        wdg = create_select(name='shape', options=options, 
                            label='', value=self.select,
                            width_label=0, width_drop=150)
        self._widgets['shape'] = wdg
        self._widgets['shape'].observe(self._on_change_shape, type='change')
        self._widgets['param'] = {}
        self._reset_widgets()

    def _create_ui(self, **kwargs):
        box = [self._widgets['shape'], self._out]
        return ipy.VBox(box)

    def _change_select(self, select):
        """Update parameters when shape type is changed."""
        self.select = select
        param = self.shapes[self.select]['param']
        self._values['param'].clear()
        self._values['param'].update(param)
        
    def _reset_widgets(self):
        """Re-fill self._widgets with input widgets for selected shape."""
        # param = self.shapes[self.select]['param']
        param = self._values['param']
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
        self._widgets['param'].clear()
        self._widgets['param'].update(items)
        self.update_values()

    def notify(self):
        """Raise notification after change.
        
        Calls `self._call_change(shape=select)`.
        """
        self._call_change(**{'shape':self.select})

    def _on_change_param(self, inp):
        key = inp['name']
        if key in self._values['param']:
            self._values['param'][key] = inp['value']
            self._call_change(**{key:inp['value']})
        else:
            raise Exception('Unknown input key: {}'.format(key))

    def _on_change_shape(self, change):
        """Process change of selected shape.
        
        Create new list of parameter widgets and display it in the output area.
        
        """
        layout = ipy.Layout(width='100%')
        if (change['name']=='value'):
            # get new set of params for the new shape
            self._change_select(change['new'])
            # create widgets for the new shape params
            self._reset_widgets()
            # display the new param widgets
            wdg = []
            for w in self._widgets['param'].values():
                wdg.append(w.ui())
            box = ipy.VBox(wdg, layout=layout)
            self._out.clear_output()
            with self._out:
                display(box)
            # callback to the UI._change method. 
            self._call_change(**{'shape':self.select})
 
    def update_widgets(self):
        """Update parameter input widgets from data."""
        if self._values['shape'] != self.select:
            self._on_change_shape({'name':'value', 'new':self._values['shape']})
        param = self._values['param']
        items = self._widgets['param']
        for key in param:
            item = items[key]
            if isinstance(item, SelectInput):
                items[key].index = param[key]
            else:
                items[key].value = param[key]          
        
    def update_values(self):
        """Update input data from widgets.
        
        Input data are stored in self._values as dict.
        
        """
        # no update if param widgets are not yet defined
        if 'param' in self._widgets and len(self._widgets['param'].keys())>0:
            self._values.clear()
            self._values['shape']=self.select
            param = {}
            items = self._widgets['param']
            for key in items:
                item = items[key]
                if isinstance(item, SelectInput):
                    param[key] = item.index
                else:
                    param[key] = item.value
            self._values['param'] = param
    
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg)
        self._on_change_shape({'name':'value', 'new':self.select})
    
    def create_shape(self):
        """Create selected shape instance from stressfit.shapes."""
        self.update_values()
        comp = shapes.create(self._values['shape'],**self._values['param'])
        return comp

        
class UI_geometry(UI_base_list):
    """Dialog for sample geometry.
    
    This UI allows to set the sample and scan geometry and mainain
    a named list of such oreiantations.
    
    Each geometry is defined by four vectors:    
    
    angles : array_like(3)
        Euler YXY angles of sample geometry
    rotctr : array_like(3)
        Rotation centre in local sample coordinates [mm]
    scandir : array_like(3)
        Scan direction in local sample coordinates
    scanorig : array_like(3)
        Scan origin (position corresponding to zero scan position)
        
    The values are returned by the method :meth:`get_values`. It returns
    named list with the keys:
        input:
            Values of the above orientation vectors from the input widgets
        list:
            A named list of oreintations added by the user. 
               
    """
    
    def __init__(self, name):
        super().__init__(name, 
                         ['angles', 'scandir', 'rotctr', 'scanorig'],
                         list_template='80px repeat(4,auto) ',
                         list_hdr=['ID','angles', 'scandir', 'rotctr', 'scanorig'],
                         list_fmt=['{}'] + 4*['({:g},{:g},{:g})'])           
        self.uiparam['title'] = 'Sample orientation'
        self.uiparam['list_title'] = 'Defined geometries'
        self.uiparam['add_label'] = 'Enter unique name (e.g. "radial", "axial", ...)'
        self.uiparam['add_label_width'] = 'auto'
        
    def _init_values(self, **kwargs):
        #print('UI_geometry._init_values\n{}'.format(kwargs))
        super()._init_values()
        # _data points to _values['list']
        self._data = self._values['list']
        self._values['input']['angles'] = [135,0,0]
        self._values['input']['rotctr'] = [0,0,0]       
        self._values['input']['scandir'] = [0,0,-1]
        self._values['input']['scanorig'] = [0,0,0]   
    
    def _create_widgets(self, **kwargs):
        inp = {}
        inp['angles'] = {'hint':'Sample rotation: YXY Euler angles [deg]'}
        inp['rotctr'] = {'hint':'Rotation centre [mm]'}       
        inp['scandir'] = {'hint':'Scan direction vector'}
        inp['scanorig'] = {'hint':'Scan origin [mm]'}        
        for k in self._keys:
            inp[k]['value'] = self._values['input'][k]
            self._widgets[k] = ArrayInput(name=k, label=k, **inp[k])

    def _create_ui(self, **kwargs):
        lst = []
        vals = list(self._widgets.values())
        for i in range(len(vals)):
            lst.append(vals[i].ui())
        return ipy.VBox(lst)

    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given ID key.
        
        By default, return a list with key followed by all items from 
        `_values['list'][key]`. Subclasses may need to modify this behaviour.
        
        """
        d = self._values['list'][key]
        rec = [key, d['angles'], d['scandir'],d['rotctr'], d['scanorig']]
        return rec
    
#    def _dbg_observe(self, inp):
#        msg = '{}:\n'.format('input')
#        msg += '{}: {}\n'.format(inp['name'], inp['value'])
#        self.message(msg, clear=True) 
#        self._prn_input(label='values', clear=False)
#        self._prn_list(label='list', clear=False)
        
#    def _prn_input(self, label='', clear=True):
#        msg = '{}:\n'.format(label)
#        msg += '{}: {}\n'.format('angles', self._values['input']['angles'])
#        self.message(msg, clear=clear) 
    
#    def _prn_list(self, label='', clear=True):
#        msg = '{}:\n'.format(label)
#        for key in self._values['list']:
#            msg += '{}: {}\n'.format(key, self._values['list'][key])
#        self.message(msg, clear=clear)     
           

class UI_sampling(UI_base_list):
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
            corresponding Sampling objects, use the method
            :meth:`get_data` for a single sampling data set.
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
    def __init__(self, name, file='', path='', nev=3000):
        hdr = ['ID','file','nrec','&lambda;','2theta','dhkl', 
               'width_x','width_y','width_z']
        fmt = ['{}', '{}', '{:d}'] + 3*['{:.5g}'] + 3*['{:.2f}']
        super().__init__(name, ['file', 'nev'],
                         list_template='80px auto repeat(7,70px)',
                         list_hdr=hdr,
                         list_fmt=fmt,
                         file=file, path=path, nev=nev) 
        
        self.uiparam['title'] = 'Sampling points'
        self.uiparam['list_title'] = 'Loaded sampling data'
        self.uiparam['add_label_width'] = 'auto'
        self.uiparam['add_button_label'] = 'Load'
        self.uiparam['add_button_icon'] = 'file'
        self.uiparam['add_button_width'] = '80px'        

    def _init_values(self, **kwargs):
        #print('UI_sampling._init_values\n{}'.format(kwargs))
        super()._init_values()
        
        f = {'file':kwargs['file'],'path': kwargs['path']}
        self._values['input']['file'] = f
        self._values['input']['nev'] = kwargs['nev']        

    def _create_widgets(self, file='', path='', nev=3000):
        # file input
        fi = FileInput(name='file', 
                       file=self._values['input']['file']['file'], 
                       path=self._values['input']['file']['path'], 
                       label='File name',
                       tooltip='File with simulated sampling points.',
                       width_label=130)
        self._widgets['file'] = fi  
        inp = ArrayInput(name='nev',
                         label='Maximum events',
                         value=self._values['input']['nev'],
                         isInt=True, step=1000, lmin=1000,
                         width_label=130)
        self._widgets['nev'] = inp        
    
    def _create_ui(self, **kwargs):
        box = [self._widgets['file'].ui(), self._widgets['nev'].ui()]
        return ipy.VBox(box)
    
    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given key."""
        s = self._data[key]
        w = s.src['width']
        rec = [key, s.file, s.src['nrec'], s.src['wav'], s.src['tth'],
                   s.src['dmean'], w[0], w[1], w[2]]
        return rec          
       
    def add_data(self, name, data=None, update_table=False):
        """Load new sampling and add it to the list."""    
        try:
            assert self._is_unique(name)
            super().add_data(name, data=data, update_table=False)
            vals = copy.deepcopy(self._values['list'][name])
            #print('UI_sampling.add_data: {}'.format(vals))
            sampling = comm.load_sampling(**vals)
            s = sampling.src
            vals['file'] = {key: s[key] for key in ['file','path']}
            vals['nev'] = s['nrec']
            self._values['input'].update(vals)
            self._values['list'][name] = vals
            self._data[name] = sampling
        except Exception as e:
            self._delete(name)
            raise(e)
        if update_table:
            self._update_table()        
        

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
    def __init__(self, name, workspace, **kwargs):
        self.wk = workspace
        super().__init__(name, ['type','value','table'],**kwargs)
        self.uiparam['title'] = 'Beam attenuation'
        self._current_att = {'att':1.1}
         
    def _init_values(self, **kwargs):
        self._values['table'] = 'Fe_mu.dat'
        self._values['value'] = 1.1
        self._values['type'] = UI_attenuation._att_types[1]
        for key in self._keys:
            if key in kwargs:
                self._values[key] = kwargs[key]

                
    def _create_widgets(self, table='Fe_mu.dat', value=1.1):        
        radio = SRadioButtons(name='type', 
                             options=UI_attenuation._att_types,
                             value=self._values['type'],
                             description='')    
        
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

    def _create_ui(self, **kwargs):
        box = []
        box.append(self._widgets['type'])
        box.append(self._widgets['value'].ui())
        box.append(self._widgets['table'].ui(width='80%'))
        return ipy.VBox(box)
    
    def _type_change(self, change):
        if (change['name']=='value'):
            t = change['new']
            is_value = (t=='value')
            self._widgets['table'].disabled = is_value
            self._widgets['value'].disabled = not is_value
                
    def get_attenuation(self):
        """Return attenuation input.
        
        Returns
        -------
        float or str
            Either single value [1/cm] or the file name.
        """
        self.update_values()
        #print('get_attenuation')
        if self._values['type'] == 'value':
            att = self._values['value']
            self._current_att.clear()
            self._current_att['att'] = att
        else:
            val = self._widgets['table'].value
            file = val['file']
            path = val['path']
            # load file only if changed
            if 'file' in self._current_att:
                qry = [val[k]==self._current_att[k] for k in ['file','path']]
            else:
                qry = [False]
            #print('get_attenuation ',qry, val)
            if not all(qry):
                att = dataio.load_data(file, kind='tables', path=path, 
                                       verbose=False)
                self._current_att.clear()
                self._current_att.update(val)
                self._current_att['att'] = att
            else:
                att = self._current_att['att']
            
        return att

    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg, **kwargs)
        self._type_change({'name':'value', 'new':self._values['type']})
        

class UI_plot_scene(UI_base):
    """Plot scene with sample geometry and sampling events.
    
    A dialog for plotting of the sample contours with scattering geometry 
    and sampling events in several projections. 
    
    Parameters
    ----------
    samplings : dict
        List of loaded sampling data sets.
        as defined by the class :class:`UI_sampling` (use :meth:`get_values`)
    geometries : dict
        List of orientations, as defined by the class :class:`UI_orientations`
        (use :meth:`get_values`)
    options : dict
        Program options
    rang : int
        Display range in mm. 
    proj : int
        Initially selected projection (0 ..2) 
    nev : int
        Number of sampling events to show.
        
    Example
    -------
    Provided that the objects `s`,`g` and 'o' are the dialogs for sampling
    , geometry and options, respectively, than create and display this dialog
    as 
    
    `plot = UI_plot_scene('plot', s.get_data(), g.get_data(), o.get_values())`
    
    `plot.show()`

    """
    
    def __init__(self, name, samplings, geometries, options, 
                 rang=16, proj=1, nev=3000):
        super().__init__(name, ['sampling','ori','nrec','proj','rang'],
                         rang=rang, proj=proj, nev=nev)    
        self.pgm_options = options

        # output area
        layout = ipy.Layout(width='75%', border='none', 
                            margin='0px 0px 0px 20px')
        self._out = ipy.Output(layout=layout) 
        # options linked to other dialogs
        self._options['sampling'] = samplings
        self._options['ori'] = geometries
        # title
        self.uiparam['title'] = 'Plot scene'
         
    def _init_values(self, **kwargs):
        self._values['nrec'] = 3000
        self._values['proj'] = 1
        self._values['rang'] = 16
        for key in self._keys:
            if key in kwargs:
                self._values[key] = kwargs[key]
                
    def _create_widgets(self, rang=16, proj=1, nev=3000):         
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
        self._widgets['rang'] = wdg 
        
    def _create_ui(self, **kwargs):
        # set selection options
        for key in self._options:
            # calling on itself makes sense: it will update widget state
            self.update_options(key, self._options[key])
        
        btn_replot = ipy.Button(description='Replot',
                                layout=ipy.Layout(width='80px'))
        btn_replot.on_click(self._on_replot)
        self._buttons.append(btn_replot)
        box1 = ipy.VBox([self._widgets['sampling'], 
                         self._widgets['ori'], 
                         self._widgets['proj'],
                         self._widgets['nrec']])
        layout = ipy.Layout(justify_content='center')
        hbox = ipy.HBox([box1, self._out],layout=layout)
        box = ipy.VBox([self._widgets['rang'], hbox])
        return box
    
    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        pfx = self.pgm_options['prefix']
        ori = self._widgets['ori'].value
        spl = self._widgets['sampling'].value
        
        if pfx:
            base = '{}_{}_{}'.format(pfx, spl, ori)
        else:
            base = '{}_{}'.format(spl, ori)
        if ext:
            fname = base + ext
        else:
            fname = base
        return fname
    
    def _on_replot(self,b):
        # Plot experiment geometry 
        # (red arrow shows motion of sampling points in stationary sample)
        self.message('')
        self.error('')
        # read selected geometry
        self.update_values()
        try:
            geometry = self._options['ori'][self._values['ori']]
            sampling = self._options['sampling'][self._values['sampling']]
        except:
            self._out.clear_output()
            self.error('Orientation or sampling is not defined.')
            return
        self._call_init(sampling=sampling, geometry=geometry)

       #  set_input(geometry, sampling)
        self._out.clear_output()
        save = self.pgm_options['save']
        if save:
            fname = self._get_output_filename(ext='.png')
        else:
            fname = ''
        with self._out:
            comm.plot_scene(self._values['nrec'], 
                            filename=fname, 
                            rang=2*[self._values['rang']],
                            proj=self._values['proj'])
            
    
class UI_resolution(UI_base):
    """Calculate resolution effects.
        
    A dialog for calculation of resolution effects: pseudo-strain, spatial
    resolution, sampling centre of mass, etc.
    
    Parameters
    ----------
    samplings : dict
        List of loaded sampling data sets.
        as defined by the class :class:`UI_sampling` (use :meth:`get_values`)
    geometries : dict
        List of oriantations, as defined by the class :class:`UI_orientations`
        (use :meth:`get_values`)
    options : dict
        Program options
    rang : int
        Display scan range in mm. 
    steps : int
        number of steps for the scan
    nev : int
        Number of sampling events to show.
        
    Example
    -------
    Provided that the objects `s`,`g` and 'o' are the dialogs for sampling
    , geometry and options, respectively, than create and display this dialog
    as 
    
    `resol = UI_resolution('resolution', s.get_data(), g.get_data(), o.get_values())`
    
    `resol.show()`
    
    
    """
    
    def __init__(self, name, samplings, geometries, options, 
                 rang=[-10, 10], steps=21, nev=3000):
        super().__init__(name, 
                         ['sampling','ori','nrec','strain','resolution',
                          'rang','steps'],
                         rang=rang, steps=steps, nev=nev) 
        self.pgm_options = options        
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none'))         
        # select options
        self._options['sampling'] = samplings
        self._options['ori'] = geometries
        # title
        self.uiparam['title'] = 'Calculate resolution effects'
         
    def _init_values(self, **kwargs):
        self._values['sampling'] = ''
        self._values['ori'] = ''
        self._values['nrec'] = kwargs['nev']
        self._values['strain'] = True
        self._values['resolution'] = True
        self._values['rang'] = kwargs['rang']
        self._values['steps'] = kwargs['steps']

    def _create_widgets(self, **kwargs): 

        # sampling selection
        wdg = create_select(name='sampling', label='Sampling', 
                            options=[], 
                            width_label=80, width_drop=100)
        self._widgets['sampling'] = wdg
        # oreiantation selection
        wdg =  create_select(name='ori', label='Orientation', 
                             options=[], 
                             width_label=80, width_drop=100)
        self._widgets['ori'] = wdg
        # number of events to plot
        wdg = create_input_int(name='nrec', label='Events',
                               value=self._values['nrec'], 
                               lmin=1000, lmax=100000, step=1000,
                               width_label=80, width_num=100)
        self._widgets['nrec'] = wdg 
        # calculate pseudo-strain
        wdg = create_checkbox(name='strain', label='Pseudo-strain', 
                              value=self._values['strain'],
                              width_label=100, width_box=20,
                              tooltip='Calculate psudo-strain',
                              indent=False)
        self._widgets['strain'] = wdg
        # calculate resolution
        wdg = create_checkbox(name='resolution', label='Resolution', 
                              value=self._values['resolution'],
                              width_label=100, width_box=20,
                              tooltip='Calculate resolution',
                              indent=False)
        self._widgets['resolution'] = wdg
        # scan range
        wdg = ArrayInput(name='rang', 
                   label='Scan range',
                   value=self._values['rang'], 
                   hint='',
                   step=0.1,
                   width_label=80)
        self._widgets['rang'] = wdg
        # number of steps
        wdg = ArrayInput(name='steps', 
                   label='Steps',
                   value=self._values['steps'], 
                   hint='',
                   isInt=True,
                   lmin=3,
                   lmax=201,
                   width_num=80,
                   width_label=80)
        self._widgets['steps'] = wdg        

    def _create_ui(self, **kwargs):
        # set selection options
        for key in self._options:
            # calling on itself makes sense: it will update widget state
            self.update_options(key, self._options[key])
        
        btn_run = ipy.Button(description='Run',
                                layout=ipy.Layout(width='80px'))
        btn_run.on_click(self._on_replot)
        self._buttons.append(btn_run)
        
        layout = ipy.Layout(margin='10px 5px 5px 10px')
        box1 = ipy.VBox([self._widgets['sampling'], 
                         self._widgets['ori'], 
                         self._widgets['nrec']],
                         layout=layout)               
        
        box2 = ipy.VBox([self._widgets['strain'], 
                         self._widgets['resolution']],
                         layout=layout)
        layout=ipy.Layout(margin='10px 5px 5px 10px', min_width='30%')
        box3 = ipy.VBox([self._widgets['rang'].ui(), 
                         self._widgets['steps'].ui()],
                         layout=layout)
        
        #style = {'margin': '5px 5px 5px 50px'}
        hbox = ipy.HBox([box1, box2, box3])
        box = ipy.VBox([hbox, self._out])
        return box
            
    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        pfx = self.pgm_options['prefix']
        ori = self._widgets['ori'].value
        spl = self._widgets['sampling'].value
        
        if pfx:
            base = '{}_{}_{}'.format(pfx, spl, ori)
        else:
            base = '{}_{}'.format(spl, ori)
        if ext:
            fname = base + ext
        else:
            fname = base
        return fname
    
    def _on_replot(self,b):
        self.message('')
        self.error('')
        # Calculate and plot results
        self.update_values()
        try:
            # read selected geometry
            geometry = self._options['ori'][self._values['ori']]
            # read selected sampling
            sampling = self._options['sampling'][self._values['sampling']]
        except:
            self._out.clear_output()
            self.error('Orientation or sampling is not defined.')
            return
        # call init function if defined
        self._call_init(sampling=sampling, geometry=geometry)

        nev = self._values['nrec']
        fname = self._get_output_filename()
        rang = list(self._values['rang'])
        nstp = self._values['steps']
        scan_range = rang + [nstp]
        save = self.pgm_options['save']
        
        # clear output
        self._out.clear_output()
        with self._out:
            if self._values['strain']:
                comm.report_pseudo_strains(scan_range, fname, 
                                           nev=nev,
                                           intensity=True,
                                           inline=True,
                                           plot=True, 
                                           save=save)
            if self._values['resolution']:    
                comm.report_resolution(scan_range, fname, 
                                       nev=nev,
                                       cog=True,
                                       inline=True,
                                       plot=True, 
                                       save=save)
      

class UI_data(UI_base_list):
    """List of input experimental data.
    
    Allows to load data sets and associate them with 
    experimental geometry and sampling events.
    
    
    The values are returned by the method :meth:`get_values`. It returns
    named list of dict objects with keys:
        file:
            File name of the data,
        data:
            Content of the datafile as numpy array
        geometry:
            Associated geometry name
        sampling:
            Associated sampling name
    
    
    Parameters
    ----------
    name : str
        Name of this UI unit.
    samplings : dict
        List of loaded sampling data sets.
        as defined by the class :class:`UI_sampling` (use :meth:`get_values`)
    geometries : dict
        List of orientations, as defined by the class :class:`UI_orientations`
        (use :meth:`get_values`)
    """

    def __init__(self, name, workspace, samplings, geometries, options, header='Input data'):
        self.wk = workspace
        super().__init__(name,  
                         ['strain', 'intensity', 'ori', 'sampling','nrec'],
                         list_template='80px auto auto repeat(2,100px)',
                         list_hdr=['ID','strain', 'intensity',  
                                   'orientation','sampling'],
                         list_fmt=['{}'] + 2*['{}'] + 2*['{}']) 
        self.pgm_options = options   
        # output area for plots
        self._outgr = ipy.Output(layout=ipy.Layout(width='100%', border='none'))
        self.uiparam['title'] = header
        self.uiparam['list_title'] = 'Loaded data'
        self.uiparam['add_label_width'] = 'auto'
        self.uiparam['add_button_label'] = 'Load'
        self.uiparam['add_button_icon'] = 'file'
        self.uiparam['add_button_width'] = '80px'        
              
        # select options
        self._options['sampling'] = samplings
        self._options['ori'] = geometries
        # title
        self.uiparam['title'] = header

    def _init_values(self, **kwargs):
        super()._init_values()
        self._values['input']['strain'] = ''
        self._values['input']['intensity'] = ''
        self._values['input']['ori'] = ''
        self._values['input']['sampling'] = ''
        self._values['nrec'] = 3000
                
    def _create_widgets(self, **kwargs):
        """Create widgets registered in self._widgets."""
        # file input
        fi = FileInput(name='strain', 
                       file='eps_SS_rad.dat', 
                       path=self.wk.full_path('data', as_posix=True),
                       fileonly=True,
                       label='Strain',
                       tooltip='File with strain data.',
                       width_label=80)
        self._widgets['strain'] = fi
        fi = FileInput(name='intensity', 
                       file='int_SS_rad.dat', 
                       path=self.wk.full_path('data', as_posix=True),
                       fileonly=True,
                       label='Intensity',
                       tooltip='File with intensity data.',
                       width_label=80)
        self._widgets['intensity'] = fi
        """
        # sampling selection
        wdg = create_select(name='sampling', label='Sampling', 
                            options=[], 
                            width_label=80, width_drop=100)
        self._widgets['sampling'] = wdg
        # geometry selection
        wdg =  create_select(name='ori', label='Orientation', 
                             options=[], 
                             width_label=80, width_drop=100)
        self._widgets['ori'] = wdg  
        """
        # sampling selection
        wdg =  SelectInput(name='sampling', label='Sampling', width_label=80, 
                           width_drop=100)
        self._widgets['sampling'] = wdg
        # geometry selection
        wdg =  SelectInput(name='ori', label='Orientation', width_label=80, 
                           width_drop=100)
        self._widgets['ori'] = wdg
        
        # number of events to plot
        wdg = create_input_int(name='nrec', label='Events',
                               value=self._values['nrec'], 
                               lmin=1000, lmax=100000, step=1000,
                               width_label=80, width_num=100)
        self._widgets['nrec'] = wdg 

    def _create_ui(self, **kwargs):
        # set selection options
        for key in self._options:
            # calling on itself makes sense: it will update widget state
            self.update_options(key, self._options[key])
        
        btn_run = ipy.Button(description='Show',
                                layout=ipy.Layout(width='80px'))
        btn_run.on_click(self._on_replot)
        self._buttons.append(btn_run)
        # dirty trick: this is not a button, but we want it next to the button
        self._buttons.append(self._widgets['nrec'])
        layout = ipy.Layout(margin='10px 5px 5px 10px')
        box = ipy.VBox([self._widgets['strain'].ui(),
                        self._widgets['intensity'].ui(),
                        self._widgets['sampling'].ui(),
                        self._widgets['ori'].ui()],
                         layout=layout)               
        return box


    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given ID key."""
        if key in self._data:   
            lst = self._values['list'][key]
            dt = self._data[key]
            ori = lst['ori']
            sam = lst['sampling']
            if dt['eps'] is not None:
                epsfile = lst['strain']
            else:
                epsfile='none'
            if dt['int'] is not None:
                intfile = lst['intensity']
            else:
                intfile='none'
        else:    
            lst = self._values['list'][key]
            ori = lst['ori']
            sam = lst['sampling']
            epsfile = 'none'
            intfile = 'none'  
        rec = [key, epsfile, intfile, ori, sam ]
        return rec 
            
    def validate(self, delete=False):
        """Check consistency of loaded data with their attributres.
        
        For example, the path or geometry names may have changed, so
        they are not valid any more. This method prints error messages for 
        all inconsistencies found. 
        
        If `delete` is True, inconsistent files are deleted 
        from the list. 
        
        """
        self.error('')
        self.update_values()
        for key in self._values['list']:
            item = self._values['list'][key]
            ans1 = item['ori'] in self._options['ori']
            if not ans1:
                msg = 'Orientation "{}" for data "{}" is not defined.'
                self.error(msg.format(item['ori'], key), clear=False)
            ans2 = item['sampling'] in self._options['sampling']
            if not ans2:
                msg = 'Sampling "{}" for data "{}" is not defined.'
                self.error(msg.format(item['sampling'], key), clear=False)
            if delete and (not (ans1 and ans2)):
                self._delete(key)

    def _get_output_filename(self, name='data_name', ext=''):
        """Generate base output filename."""
        pfx = self.pgm_options['prefix']        
        if pfx:
            base = '{}_{}'.format(pfx, name)
        else:
            base = '{}'.format(name)
        if ext:
            fname = base + '_{}'.format(ext)
        else:
            fname = base
        return fname        
        
    def add_data(self, name, data=None, update_table=False):
        """Add a dataset to the list.
        
        At the end, calls self._call_change with the list content.
        
        Parameters
        ----------
        name : str
            ID string for new item. It must be a unique ID string.
        data : dict
            New keys with values to be added to the table.
            If None, then values are taken from corresponding widgets.
        update_table : bool
            If True, call self._update_table() to redraw the table content.
        """

        try:
            assert self._is_unique(name)
            self._dbg_message('UI_data.add_data data: {}'.format(data))

            super().add_data(name, data=data, update_table=False)
            vals = copy.deepcopy(self._values['list'][name])
            self._dbg_message('UI_data.add_data values: {}'.format(vals))
            ori = vals['ori']
            sam = vals['sampling']
            ifile = vals['intensity'].strip()
            sfile = vals['strain'].strip()  
            #print('add_data '+ifile)
            if not ifile or ifile=='.':
                ifile = None
            try:
                geom = self._options['ori'][ori]
                sampling = self._options['sampling'][sam]
            except:
                self._outgr.clear_output()
                self.error('Cannot associate orientation or sampling with the data.')
                return
            content = comm.load_input(sfile, intensity=ifile,
                                      path=None, 
                                      angles=geom['angles'],
                                      rotctr=geom['rotctr'],
                                      scandir=geom['scandir'],
                                      scanorig=geom['scanorig'],
                                      sampling=sampling, 
                                      verbose=False)
            self._dbg_message('UI_data.add_data content: {}'.format(content))
            self._data[name] = content
            if update_table:
                self._update_table()
        except Exception as e:
            self._delete(name)
            raise(e)
    
    def _plot_comparison(self, simdata, expdata):
        
        def exp_to_dict(expdata, what='int'):
            # convert exp. data do dict for plotting
            out = {}
            if what=='int':
                out['title'] = expdata['intfile']
                out['xlabel'] = 'Scan position, mm'
                out['ylabel'] = 'Intensity, rel. units'
                out['x'] = expdata['int'][:,0]
                out['y'] = expdata['int'][:,1]
                out['yerr'] = expdata['int'][:,2]
            else:
                out['title'] = expdata['epsfile']
                out['xlabel'] = 'Scan position, mm'
                out['ylabel'] = 'Strain,  1e-6'
                out['x'] = expdata['eps'][:,0]
                out['y'] = expdata['eps'][:,1]
                out['yerr'] = expdata['eps'][:,2]
            return out
        
        def scale(data, A=1, B=0, x0=0):
            """Scale intensity data."""
            out = copy.deepcopy(data)
            out['x'] = data['x']-x0
            out['y'] = A*data['y'] + B
            if 'yerr' in data:
                out['yerr'] = A*data['yerr']
            return out
        
        coll = {}
        has_intensity = None
        dim = 2
        inline = True
        for key in expdata:
            dexp = expdata[key]
            dsim = simdata[key]  
            # use first data set to decide if intenisty should be plotted
            if has_intensity is None:
                has_intensity = dsim['intensity'] is not None
                if has_intensity:
                    dim = 2
                    inline = True
            # strain
            toplot = exp_to_dict(dexp, what='eps')
            toplot['title'] = dexp['epsfile']
            toplot['label'] = 'experiment'
            toplot['fmt'] = 'ko'            
            other = dsim['strain']
            other['label'] = 'simulation'
            other['fmt'] = 'b-'
            toplot['other'] = other
            coll[key+'_s'] = toplot            
            # intensity
            if dexp['int'] is not None:
                toplot = exp_to_dict(dexp, what='int')
                toplot['title'] = dexp['intfile']
                toplot['label'] = 'experiment'
                toplot['fmt'] = 'ko'
                sim_int = dsim['intensity']['y']
                A = dexp['int'][:,1].mean()/sim_int.mean()
                other = scale(dsim['intensity'], A=A, B=0, x0=0)
                other['label'] = 'simulation'
                other['fmt'] = 'b-'
                toplot['other'] = other
            else:
                toplot = dsim['intensity']
                toplot['label'] = 'simulation'
                toplot['fmt'] = 'b-'
            coll[key+'_i'] = toplot
        if len(coll)>0:    
            gr.plot_collection(coll, dim=dim, inline=inline)      
    
    def _on_replot(self,b):
        """Plot data with simulated pseudo-strains and pseudo-intensities."""
        expdata = {}
        simdata = {}
        self.message('')
        self.error('')
        for name in self._data:
            data = self._data[name]
            x = data['eps'][:,0]        
            scan_range = [min(x), max(x), 2*len(x)+1]
            fname = self._get_output_filename(name=data['epsfile'])       
            nev = self._values['nrec']
            save = self.pgm_options['save']
            
            # read selected geometry
            geometry = {k:data[k] for k in ['angles','scandir','rotctr','scanorig']}
            # read selected sampling
            sampling = data['sampling']
            # call init function if defined
            self._call_init(sampling=sampling, geometry=geometry)
            with self._msg:
                res = comm.report_pseudo_strains(scan_range, fname, 
                                                 nev=nev,
                                                 intensity=True,
                                                 inline=True,
                                                 plot=False, 
                                                 save=save)
            expdata[name] = data
            simdata[name] = res
        # clear output
        self._outgr.clear_output()
        with self._outgr:
            self._plot_comparison(simdata, expdata)
      
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg, **kwargs)
        display(self._outgr)
        
#%% Top level UI class

class UI():
    """Top level class generating the Jupyter notebook interface."""
    
    _registered_ui = ['workspace',
                      'options',
                      'shape',
                      'geometry',
                      'sampling',
                      'attenuation',
                      'scene',
                      'resolution',
                      'data']
    def __init__(self):
        self._last_input_file = 'input.json'
        self.setup = {}
        self.ui = {}
        self.wk = dataio.workspace()
        self._msg = ipy.Output(layout=ipy.Layout(width='100%'))
        self._err = ipy.Output(layout=ipy.Layout(width='100%'))
        # setup dialog
        self._add_input_ui(UI_workspace('workspace'))
        self._add_input_ui(UI_options('options'))
        self._add_input_ui(UI_shape('shape',select=shapes.Tube))
        self._add_input_ui(UI_geometry('geometry'))
        self._add_input_ui(UI_sampling('sampling'))
        self._add_input_ui(UI_attenuation('attenuation', self.wk))
        
        # initialize
        self.ui['geometry'].add_data('radial', update_table=False)
        inpdir = self.ui['workspace'].get_values()['data']
        inp = {'input':{}}
        inp['input']['file'] = {'path':inpdir, 'file':'events_S_1mm.dat'}
        inp['input']['nev'] = 3000
        self.ui['sampling'].set_values(inp)
        #print(self.ui['sampling'].get_values())
        self.ui['sampling'].add_data('1x1x5mm', data=inp['input'], update_table=False)
        
        # execution dialogs
        obj = UI_plot_scene('scene',
            self.ui['sampling'].get_data(),
            self.ui['geometry'].get_data(),
            self.ui['options'].get_values())
        self._add_input_ui(obj)
        
        obj = UI_resolution('resolution',
            self.ui['sampling'].get_data(),
            self.ui['geometry'].get_data(),
            self.ui['options'].get_values())
        self._add_input_ui(obj)
        
        obj = UI_data('data',
            self.wk,
            self.ui['sampling'].get_data(),
            self.ui['geometry'].get_data(),
            self.ui['options'].get_values())
        # initialize dat list
        inp = {'strain':'eps_SS_rad.dat', 
               'intensity':'int_SS_rad.dat',
               'ori':'radial',
               'sampling':'1x1x5mm'}
        obj.set_values({'input':inp})
        obj.add_data('radial', update_table=False)
        obj.validate()
        self._add_input_ui(obj)
        
        # set event handlers
        self.ui['workspace'].set_on_change(self._change)
        self.ui['shape'].set_on_change(self._change)
        self.ui['geometry'].set_on_change(self._change)
        self.ui['sampling'].set_on_change(self._change)
        
        self.ui['scene'].set_on_init(self._init)
        self.ui['resolution'].set_on_init(self._init)
        self.ui['data'].set_on_init(self._init)
    
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
        # Define tabs_data in the order of tabs
        layout=ipy.Layout(width='100%')        
        tabs_data = {}
        tabs_data['workspace'] = {'title':'Workspace',
                                  'ui':[self.ui['workspace'],
                                        self.ui['options']]}
        tabs_data['geometry'] = {'title':'Geometry',
                                  'ui':[self.ui['shape'],
                                        self.ui['geometry']]}
        tabs_data['material'] = {'title':'Material',
                                  'ui':[self.ui['attenuation']]}
        tabs_data['sampling'] = {'title':'Sampling',
                                  'ui':[self.ui['sampling'],
                                        self.ui['scene']]}
        tabs_data['resolution'] = {'title':'Resolution',
                                  'ui':[self.ui['resolution']]}
        tabs_data['data'] = {'title':'Data',
                                  'ui':[self.ui['data']]}
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
        
        # button bar
        btn_save = ipy.Button(description='Save input')
        btn_save.on_click(self._on_save_input)
        btn_load = ipy.Button(description='Load input')
        btn_load.on_click(self._on_load_input)
        
        display(ipy.HBox([btn_save,btn_load]))
        
        display(tab)
        keys = list(tabs_data.keys())
        for key in keys:
            for ui in tabs_data[key]['ui']:
                with tabs[key]:
                    ui.show(msg=self._msg, err=self._err)
        

        display(self._err)
        display(self._msg) 
    
    def _on_save_input(self,b):
        s = choose_file_save(initialdir=self.wk.path('work').as_posix(), 
                        initialfile=self._last_input_file,
                        filetypes=(('Setup files','*.json')))
        if s:
            self.save(filename=s)
            p = _Path(s)
            self._last_input_file = p.name
            
    def _on_load_input(self,b):
        s = choose_file(initialdir=self.wk.path('work').as_posix(), 
                        initialfile=self._last_input_file,
                        filetypes=(('Setup files','*.json')))
        if s:
            self.load(filename=s)
            p = _Path(s)
            self._last_input_file = p.name
            
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
                self.ui['scene'].update_options('ori',new_options)
                self.ui['resolution'].update_options('ori',new_options)
                self.ui['data'].update_options('ori',new_options)
                self.ui['data'].validate()
        elif obj.name in ['sampling']:
                data.update(**kwargs)
                #new_options = self.ui['sampling'].get_list()
                new_options = self.ui['sampling'].get_data()
                #print('UI._change, new_options:\n{}'.format(new_options))
                self.ui['scene'].update_options('sampling',new_options)
                self.ui['resolution'].update_options('sampling',new_options)
                self.ui['data'].update_options('sampling',new_options)
                self.ui['data'].validate()

    def _init(self, obj, **kwargs):
        """Initialize stressfit.
        
        Callback used by IU groups to initialize stressfit.
        """ 
        if obj.name in ['scene', 'resolution', 'data']:
            self._err.clear_output()
            if 'sampling' in kwargs:
                #print('UI._init: {}'.format(kwargs))
                comm.set_sampling(kwargs['sampling'])
            if 'geometry' in kwargs:
                comm.set_geometry(kwargs['geometry'])
            # set attenuation
            att = self.ui['attenuation'].get_attenuation()
            comm.set_attenuation(att)
            

    def save(self, filename=''):
        """Save input data in JSON format."""
        try:
            for key in self.ui: 
                self.ui[key].update_values()
            out = {'ui': self.setup} 
            txt = json.dumps(out,indent=4)
            if filename:
                #wkp = self.wk.get_paths(keys=['work'])
                #file = wkp['work'].joinpath(filename)
                f = open(filename, 'w')
                f.write(txt)
                f.close()
            else:
                with self._msg:
                    print(txt)
        except Exception as e:
            with self._err:
                print(e)
            raise e
        
    def load(self, filename=''):
        """Load input data in JSON format."""
        try:
            f = open(filename, 'r')
            lines=f.readlines()
            f.close()
            inp = json.loads('\n'.join(lines))
            if not 'ui' in inp:
                with self._err:
                    print('Wrong setup data format.')
                return
            for key in self.ui:
                if key in inp['ui']:
                    #print('updating {}'.format(key))
                    self.ui[key].set_values(inp['ui'][key])
                    self.ui[key].update_widgets() 
                    self.ui[key].notify()
            self.ui['workspace'].validate_workspace()
        except Exception as e:
            with self._err:
                print(e)
            raise e

