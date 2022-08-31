"""Classes and functions for creating notebook interface.

Created on Tue Oct  5 11:38:15 2021
@author: Jan Saroun, saroun@ujf.cas.cz
"""

# TODO manage messages in widgets and config, define extra spaces for it in notebook
# TODO model definition - strain
# TODO fit strain
# TODO input - compliance
# TODO model definition - stress
# TODO fit stress
# 

try:
    from pygments import formatters, highlight, lexers
    _has_pygments = True
except ImportError:
    _has_pygments = False

#print('_has_pygments = {}'.format(_has_pygments))
#import abc
import re
import traceback
import abc
import ipywidgets as ipy
import json
import copy
from IPython.display import display
from .widgets import DirInput, FileInput, ArrayInput, SelectInput, FitControl
from .widgets import create_header, create_select 
from .widgets import create_checkbox, create_text
from .widgets import SButton, SRadioButtons, DistTable, ScaleTable
from .widgets import choose_file, choose_file_save
from pathlib import Path as _Path
import numpy as np
import stressfit.commands as comm
import stressfit.shapes as shapes
import stressfit.dataio as dataio
import stressfit.graphs as gr
import stressfit.mccfit as mc
from stressfit.ui.config import uiconfig
from stressfit.geometry import Geometry


# get ui configuration and data with default input values
_uiconf = uiconfig()

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
        self._values = {} # dict of widget values links to _uiconf by assign
        self._buttons = [] # Command buttons next to the collection label 
        self._on_change = None
        self._options = {} # containes links to selection lists
        self._init_values(**kwargs)
        
# abstract methods


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

    def _clean_keys(self, data:dict):
        """Delete all data items with undefined keys."""
        for key in list(data.keys()):
            if not key in self._keys:
                del data[key]

    def _check_keys(self,data):
        """Verify that data contains all required keys."""
        if not _has_keys(data,self._keys):
            msg = 'Required key(s) missing in the input data: {}'
            raise Exception(msg.format(self._keys))
 
    def _call_change(self, **kwargs):
        """Call when the UI changes state."""
        if self._on_change:
            self._on_change(self, **kwargs)      

    def _init_values(self, **kwargs):
        """Link values to global input data.
        
        kwargs may contain parameters which replace the global setting.
        
        Called by constructor. 
        """
        self.assign()
        for key in kwargs:
            if key in self._keys:
                _uiconf.update(key,kwargs[key])
            
# properties        
    @property
    def name(self):
        """ID string for the input block."""
        return self._name


# public methods
    def assign(self):
        """Assign pointers to configuration data to self._values."""
        self._values = _uiconf.get(self.name)

    def notify(self):
        """Raise notification after change.
        
        UI_base does nothing. Subclasses may need to implement its own actions.
        """

    def set_values(self, values):
        """Set values.
        
        Write values into configuration data (_uiconf). 
        To update widgets for the new values, call :meth:`update_widgets`.
        
        Parameters
        ----------
        values : dict
            New values.
        """
        vals = copy.deepcopy(values)
        self._clean_keys(vals)
        _uiconf.update(self.name, vals)

    def get_values(self):
        """Return actual input values as dict."""
        self.update_values()
        return self._values
    
    def update_widgets(self):
        """Update input widgets from values."""        
        for key in self._widgets:
            if key in self._values:
                self._widgets[key].value = copy.deepcopy(self._values[key])

    def update_values(self):
        """Update values from widgets."""
        vals = {}
        for key in self._keys:
            if key in self._widgets:
                vals[key] = self._widgets[key].value
        self.set_values(vals)

    def update_options(self, name):
        """Update the selection list.
        
        Parameters
        ----------
        name : str
            Name of the options list
        """
        if not name in self._options:
            return
        if not name in self._widgets:
            msg = 'No selection widget of this name: {}'
            raise Exception(msg.format(name))    
        self._widgets[name].options = self._options[name]

    def message(self, txt, clear=True):
        """Print message in the message output, if defined."""
        if self._msg:
            if clear:
                self._msg.clear_output()
            with self._msg:
                print(txt)
        else:
            print(txt)
                
    def error_trace(self, clear=True):
        """Format and print traceback from the last exception."""        
        if self._err and clear:
            self._err.clear_output()
        tb_text = "".join(traceback.format_exc()) 
        if _has_pygments:
            lexer = lexers.get_lexer_by_name("pytb", stripall=True)
            if self._err:
                style = '<STYLE>\n'
                style += '.highlight .gr { color: #FF0000 }\n'
                style += '.highlight .o { color: #333333 }\n'
                style += '.highlight .m { color: #6600EE; font-weight: bold }\n'
                style += '.highlight .nb { color: #007020 }\n'
                style += '.highlight .gt { color: #0044DD }\n'
                style += '.highlight .s1 { background-color: #fff0f0 } \n'
                style += '.highlight .bp { color: #902020 } \n'
                style += '.highlight .k { color: #008800; font-weight: bold }\n'   
                style += '</STYLE>\n'
                formatter = formatters.get_formatter_by_name('html',style='colorful')
                tb_colored = highlight(tb_text, lexer, formatter)
                t = ipy.HTML(value=style+tb_colored)
                with self._err:
                    display(t)
            else:
                self.message('terminal16m')
                formatter = formatters.get_formatter_by_name('terminal16m')
                tb_colored = highlight(tb_text, lexer, formatter)
                print(tb_colored)
        elif self._err:
            with self._err:
                print(tb_text)
        else:
            self.message('no _err')
            print(tb_text)
        
    def error(self, txt, clear=True, trace=True):
        """Print error message in the error output, if defined."""
        if clear and self._err:
            self._err.clear_output()
        if not txt:
            return
        if self._err:
            with self._err:
                if isinstance(txt, Exception):
                    print(txt)
                else:
                    s1 = "<font color='red'>{}</font>".format(txt)
                    s = re.sub(r'\n',r'<br/>',s1)
                    t = ipy.HTML(value=s)
                    display(t)
        else:
            print(txt)
        if trace:
            self.error_trace(clear=False)

    def set_on_change(self, on_change):
        """Set callback to be executed when the UI changes state."""
        self._on_change = on_change
              
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
        self.uiparam['add_name_width'] = '120px'
        self.uiparam['add_label'] = 'Name'
        self.uiparam['add_label_width'] = 'auto'
        
        # save format info for the table of data sets.
        self._template = list_template
        self._hdr = list(list_hdr)
        self._fmt = list(list_fmt)
    
    def _redraw_table(self):
        """Update stressfit configuration and redraw the table of data sets."""
        try: 
            grid = self._create_list()      
            self._out.clear_output()
            if len(self._values['list'].keys())>0:
                    with self._out:
                        display(grid)                    
        except Exception:
            msg = 'Cannot redraw table.'
            self.error(msg, trace=True, clear=False)
            pass

    def _delete(self, key):
        """Delete given item from the list.
        
        If successful, updates the table and calls :meth:`notify`.
        """
        if key in self._values['list']:
            _uiconf.data.delete(self.name, key)
            self._redraw_table()
            self.notify()
            
    def notify(self):
        """Raise notification after change.
        
        Calls `self._call_change(list=self._values['list'])`.
        """
        arg = {'list':self._values['list']}
        #self.message('notify {}'.format(self.name))
        self._call_change(**arg)
        
    def _on_del_click(self, b):
        """Delete the table row using `delete` button."""
        self._delete(b.value)
            

    def _on_edit_click(self, b):
        """Update input widgets with values from selected list item."""
        data = copy.deepcopy(self._values['list'][b.value])
        self.set_values({'input':data})
        self._name_inp.value = b.value
        self.update_widgets()

    def _on_add(self,b):
        """Call :meth:`add_data`.

        `on_change` is sent at the end with 'list' and 'data' arguments.
        """
        txt = self._name_inp.value
        
        self._dbg_message('Add {}'.format(txt), clear=True)
        
        self.add_data(txt, update_ui=True)
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
        grid_layout = ipy.Layout(grid_template_columns='50px 50px ' + self._template)  
        hdr = ['',''] + self._hdr
        fmt = self._fmt
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
            be = SButton(description='',icon='edit',value=key,
                            layout=ipy.Layout(width='30px'), tooltip='Modify')
            be.on_click(self._on_edit_click)
            wrec = [b, be]
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

    def _init_values(self, **kwargs):
        """Link values to global input data.
        
        kwargs may contain parameters which replace the global setting.
        
        Called by constructor. 
        """
        self.assign()
        
    def assign(self):
        """Assign values to global data."""
        self._values['list'] = _uiconf.data.get_item(self.name)
        self._values['input'] = _uiconf.data.inputs[self.name]
        self._values['conf'] = _uiconf.get(self.name)
        self._keys_conf = []

    def set_values(self, values, item=None):
        """Set values.
        
        Set new input data. 
        
        To update also widgets for the new values, call :meth:`update_widgets`.
        
        Parameters
        ----------
        values : dict
            Input data. The expected top-level keys are:
            
            input : str, optional
                Keys and values for the input widgets.
            list : dict, optional
                Data sets to be added in the table. Each data set 
                is given a unique key in values['list']. If item is
                defined, values must contaiin only that item values.
        item : str
            Item key to be updated for data list. If None, 
            updates all items given in values.
        """
        if 'input' in values:
            vals = copy.deepcopy(values['input'])
            self._clean_keys(vals)
            _uiconf.data.inputs[self.name].update(vals) 
        elif 'list' in values: 
            vals = copy.deepcopy(values['list'])
            for k in vals:
                self._clean_keys(vals[k])
            _uiconf.data.update(self.name, vals, item=item) 
        elif 'conf' in values: 
            vals = copy.deepcopy(values['conf'])
            vals_filt = {k:vals[k] for k in self._keys_conf}
            _uiconf.update(self.name, vals_filt)
            
    def update_values(self):
        """Update input data from widgets."""
        vals = {}
        for key in self._keys:
            if key in self._widgets:
                vals[key] = self._widgets[key].value
        self.set_values({'input':vals})
        if len(self._keys_conf) > 0:
            vals = {}
            for key in self._keys_conf:
                if key in self._widgets:
                    vals[key] = self._widgets[key].value
            self.set_values({'conf':vals})
        
    def update_widgets(self):
        """Update input widgets from data."""
        #if self.name=='imodel':
        #    self.message('{} update_widgets'.format(self.name), clear=False)
        for key in self._keys:
            if key in self._widgets:
                self._widgets[key].value = copy.deepcopy(self._values['input'][key])
        self._redraw_table()
        for key in self._keys_conf:
            if key in self._widgets:
                self._widgets[key].value = copy.deepcopy(self._values['conf'][key])

    def add_data(self, item, data=None, update_ui=False):
        """Add a dataset to the list.
        
        If an item already exists, replace it with new data.
        
        By default, put a record with all widgets values on the list.
        Subclasses may need to override this e.g. by loading data files etc.
        
        Widgets are not updated unless required by the update_ui flag.  
        Note that widgets may not exist when this method is called. 
        
        Parameters
        ----------
        item : str
            ID string for the item.
        data : dict
            Data set to be added to the table.
            If None, then values are taken from corresponding widgets.
        update_ui : bool
            If True, call :meth:`update_widgets` to redraw the ui.
        """
        if item is None or not item.strip():
            self.error('Give an ID string for the list item.')
            return
        # no input data provided: get it from the widgets content
        if data is None:
            self.update_values()         
        # else set the input data as current input values 
        else:
            try:
                self._check_keys(data)
                self.set_values({'input':data})
            except Exception as e:
                self.error(e)
        # add or replace the named item on the list using input values
        new_input = {item: self._values['input']}
        
        self.set_values({'list':new_input})           
        # update the table with data list now if required
        if update_ui:
            self.update_widgets()
        
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
        add_box = ipy.HBox([btn_add, self._name_inp, hint], 
                           layout=ipy.Layout(margin='1.5ex 0ex 0ex 0ex'))
        
        # display everything
        top = [hdr] + self._buttons
        top_wdg = ipy.HBox(top)
        layout = ipy.Layout(margin='0px 0px 20px 0px', width='100%')
        display(ipy.VBox([top_wdg, ui, add_box, list_hdr, self._out], 
                         layout=layout))
        # update table
        self._redraw_table()      
        
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
        super().__init__(name, self.wk.keys, exclude=exclude)
        self.uiparam['title'] = 'Workspace'
        self._out = ipy.Output(layout=ipy.Layout(width='100%'))     
        
    def _init_values(self, **kwargs):
        if 'exclude' in kwargs:
            self._excluded_keys = kwargs['exclude'] 
        else:
            self._excluded_keys = []
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
    
    def _create_ui(self, **kwargs):    
        rst = ipy.Button(description='Reset',
                         tooltip='Reset workspace configuration to package default.',
                         layout=ipy.Layout(width='10%'))
        rst.on_click(self._on_reset)
        save = ipy.Button(description='Save',
                          tooltip='Save workspace configuration.',
                          layout=ipy.Layout(width='10%'))
        save.on_click(self._on_save)
        load = ipy.Button(description='Reload',
                          tooltip='Re-load workspace configuration.',
                          layout=ipy.Layout(width='10%'))
        load.on_click(self._on_reload)
        
        layout = ipy.Layout(display='flex',flex_flow='row',
                                border='none',width='100%')
        appl_box = ipy.HBox([save, load, rst], layout=layout)
        
        lst = []
        for w in self._widgets:
            lst.append(self._widgets[w].ui())
        lst.append(appl_box)
        lst.append(self._out)
        return ipy.VBox(lst)
        
    def set_from_wks(self):
        """Update UI with actual workspace settings."""
        self._values.update(self.wk.get_paths())
        self.update_widgets()
        self.update_values()
    
    def _on_path_change(self, inp):
        key = inp['name']
        value = inp['value']
        self._out.clear_output()
        with self._out:
            if key == 'work':
                self.wk.change_workspace(root=value, create_tree=False, 
                                         save=False, verbose=True)            
            else:
                self.wk.set_paths(**{key:value})
            self.set_from_wks()
    
    def _on_reset(self, b):
        """Set workspace configuration to package default."""
        self.reset_workspace()
        
    def _on_save(self, b):
        """Save workspace configuration in workspace root directory."""
        self._out.clear_output()
        try:
            with self._out:
                self.wk.save()
                self.wk.validate_paths()
                print('Workspace configuration saved.')
        except Exception as e:
            self.error(e)
    
    def _on_reload(self, b):
        """Re-load workspace configuration from workspace root directory."""
        inp = {'name':'work', 'value':self._values['work']}
        try:
            self._on_path_change(inp)
        except Exception as e:
            self.error(e)      
                            
    def assign(self):
        """Do nothing.
        
        Workspace is not linked to ui_config data. 
        """
        pass
        
    def update_workspace(self, validate=True, verbose=True):
        """Update workspace configuration according to _values.
        
        If validate, check workspace configuration (paths must exist).
        """
        self._out.clear_output()
        with self._out:
            try:
                self.wk.change_workspace(root=self._values['work'], 
                                         create_tree=False, 
                                         save=False, 
                                         verbose=verbose)
                self.wk.set_paths(**self._values)
                if validate:
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
            self.set_from_wks()


class UI_options(UI_base):
    """UI for setting program options."""    
    
    def __init__(self, name, prefix='', save=True):
        super().__init__(name, ['prefix', 'save'],
                         prefix=prefix, save=save)         
        self.uiparam['title'] = 'Options'         
        
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
        return ipy.HBox([self._widgets['save'], self._widgets['prefix']])
    
    def _observe(self, change):
        if change['name']=='value':
            self._values[change['owner'].name] = change['new']

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
    
    def __init__(self, name):
        self.shapes = dataio.load_config('shapes')
        super().__init__(name, ['shape', 'param'])        
        self.uiparam['title'] = 'Sample shape'
        self._out = ipy.Output(layout=ipy.Layout(border='1px solid'))
        self._shape_info = None

    def _init_values(self, **kwargs):
        self.assign()
        self.select = self._values['shape']
        self._param = self.shapes[self.select]['param']
               
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
        self._param = self.shapes[self.select]['param']
        self._values['param'].clear()
        for key in self._param:
            self._values['param'][key] = self._param[key]['value']
        
    def _on_load(self, obj):
        try:
            self.update_values()
            fn = self._values['param']['filename']
            fullname = dataio.get_input_file(fn['file'], path=fn['path'])
            obj = shapes.create(self._values['shape'],filename=fullname)
            comm.set_shape(obj)
            
            # print info about loaded shape
            if self._shape_info is not None:
                self._shape_info.clear_output()
                with self._shape_info:
                    print(type(obj).__name__)
                    print(obj.get_param())
                
        except Exception as e:
            self.error(e)
        
    def _reset_widgets_file(self):
        """Create widgets for file input, if select = File."""
        # param = self.shapes[self.select]['param']
        param = self._param # self._values['param']
        self._widgets['shape'].value = self.select
        items = {}
        for key in param:
            p = param[key]
            value =p['value'] 
            desc = p['label']
            hint = p['hint']
            if isinstance(value, str):
                wgt = FileInput(name=key, 
                       file=value, 
                       path=dataio.workspace().path('data'),
                       #fileonly=True,
                       filetypes=(('Setup files','*.json')),
                       label=desc,
                       tooltip=hint,
                       width='50%',
                       width_label=130)
            else:
                raise Exception('Unknown parameter type: {}'.format(value))
            wgt.observe(self._on_change_param)
            items[key] = wgt
        self._widgets['param'].clear()
        self._widgets['param'].update(items)
        self.update_values()
        
    def _reset_widgets(self):
        """Re-fill self._widgets with input widgets for selected shape."""
        # param = self.shapes[self.select]['param']
        if self.select == 'File':
            self._reset_widgets_file()
            return
        param = self._param # self._values['param']
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
                wgt = ArrayInput(name=key, value=value, label=desc, hint=hint,
                                 isInt = isinstance(value, int))
            elif isinstance(value, list):
                if isinstance(value[0],float):
                    wgt = ArrayInput(name=key, value=value, label=desc, hint=hint)
                else:
                    wgt = SelectInput(name=key, options=value, 
                                      index=0,
                                      label=desc, 
                                      hint=hint,
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
        if self.select == 'File':
            self._on_load(None)

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
                
            if self.select == 'File':
                appl = ipy.Button(description='Load',
                          tooltip='Load shape from the file.',
                          layout=ipy.Layout(width='10%'))
                appl.on_click(self._on_load)
                wdg.append(appl)
                self._shape_info = ipy.Output(layout=ipy.Layout(border='none'))
                wdg.append(self._shape_info)
            else:
                self._shape_info = None
                
            box = ipy.VBox(wdg, layout=layout)
            self._out.clear_output()
            with self._out:
                display(box)
            # callback to the UI._change method. 
            self._call_change(**{'shape':self.select})
 
    def assign(self):
        """Assign values to global data."""
        self._values = _uiconf.data.get_item(self.name)   
 
    def update_widgets(self):
        """Update parameter input widgets from data."""
        if self._values['shape'] != self.select:
            # NOTE: change_shape resets parameers to default
            # we have to restore them to current values
            param = copy.deepcopy(self._values['param'])
            self._on_change_shape({'name':'value', 'new':self._values['shape']})
            self._values['param'].clear()
            self._values['param'].update(param)
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
        # trigger displaying of shape parameters
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
                         list_template='120px repeat(4,auto) ',
                         list_hdr=['name','angles', 'scandir', 'rotctr', 'scanorig'],
                         list_fmt=['{}'] + 4*['({:g},{:g},{:g})'])           
        self.uiparam['title'] = 'Sample orientation'
        self.uiparam['list_title'] = 'Defined geometries'
        self.uiparam['add_label'] = 'Enter unique name (e.g. "radial", "axial", ...)'
            
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

class UI_sampling(UI_base_list):
    """Input for a list of simulated sampling data sets.
    
    This UI permits to define a file with simulated sampling events and
    load required number of events from it. In addition, this UI 
    maintains a list of already loaded sampling files. 
    
    The values are returned by the method :meth:`get_values`. It returns
    named list with the keys:
        input:
            Filename, path and number of events for loading the sampling events . 
        list:
            A named list of sampling data sets added by the user.
            Only meta-data are included (file, path, nrec).
    
    Parameters
    ----------
    file : str
        Initial sampling file name
    path : str
        Initial directory name with sampling files
    nrec : int
        Initial number of sampling points to be loaded.
    """
    
    style_lbl = {'description_width': '100px'}
    def __init__(self, name, file='', path='', nrec=3000):
        self.wk = dataio.workspace()
        hdr = ['name','file','nrec','&lambda;','2theta','dhkl','width_xyz']
        fmt = ['{}', '{}', '{:d}'] + 3*['{:.5g}'] + ['({:.2f}, {:.2f}, {:.2f})']
        super().__init__(name, ['file', 'nrec'],
                         #list_template='80px auto repeat(7,70px)',
                         list_template='80px auto repeat(4,70px) auto',
                         list_hdr=hdr,
                         list_fmt=fmt,
                         file=file, path=path, nrec=nrec) 
        
        self.uiparam['title'] = 'Sampling points'
        self.uiparam['list_title'] = 'Loaded sampling data'
        self.uiparam['add_button_label'] = 'Load'
        self.uiparam['add_button_icon'] = 'file'
        self.uiparam['add_button_width'] = '80px'
        self.uiparam['add_label'] = 'Sampling name'

    def _create_widgets(self, file='', path='', nrec=3000):
        # file input
        fi = FileInput(name='file', 
                       file=self._values['input']['file']['file'], 
                       path=self.wk.full_path('data', as_posix=True), 
                       label='File name',
                       tooltip='File with simulated sampling points.',
                       width_label=130)
        self._widgets['file'] = fi  
        inp = ArrayInput(name='nrec', label='Maximum events', hint='',
                         value=self._values['input']['nrec'],
                         isInt=True, step=1000, lmin=1000,
                         width_label=130)
        self._widgets['nrec'] = inp        
    
    def _create_ui(self, **kwargs):
        box = [self._widgets['file'].ui(), self._widgets['nrec'].ui()]
        return ipy.VBox(box)
    
    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given key."""
        s = _uiconf.data.get_item('sampling', item=key)['fdata']
        try:
            rec = [key, s.file, s.src['nrec'], s.src['wav'], s.src['tth'],
                       s.src['dmean'], list(s.src['width'])]
        except:
            rec = [key, 'empty'] + [0] + 3*[0.0] + [3*[0.0]]
        return rec          
      
        
    def add_data(self, item, data=None, update_ui=False):
        """Load new sampling and add it to the list."""    
        try:
            super().add_data(item, data=data, update_ui=False)
            vals = copy.deepcopy(self._values['list'][item])
            #self._values['input'].update(vals)
            _uiconf.data.load(self.name, item=item)
            sampling = _uiconf.data.get_item(self.name, item=item)['fdata']
            s = sampling.src
            vals['file'] = {key: s[key] for key in ['file','path']}
            vals['nrec'] = s['nrec']
            self._values['input'].update(vals)
            self._values['list'][item].update(vals)           
        except Exception:
            if item in self._values['list']:
                self._delete(item)
            self.error_trace(clear=False)
        if update_ui:
            self.update_widgets()        
        

class UI_attenuation(UI_base):
    """Dialog for beam attenuation.

    Parameters
    ----------
    data: dict
        Initial input values (optional)
    
    """
    
    _att_types = ['value', 'table']
    def __init__(self, name, **kwargs):
        self.wk = dataio.workspace()
        super().__init__(name, ['type','value','table'],**kwargs)
        self.uiparam['title'] = 'Beam attenuation'
        self._current_att = {'att':1.1}
                
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
                       file=self._values['table']['file'], 
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

    def assign(self):
        """Assign values to global data."""
        self._values = _uiconf.data.get_item(self.name)
        
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg, **kwargs)
        self._type_change({'name':'value', 'new':self._values['type']})
        

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
    header : str
        Title of the input panel.

    """

    def __init__(self, name, header='Input data'):
        self.wk = dataio.workspace()
        super().__init__(name,  
                         ['strain', 'intensity', 'geometry', 'sampling'],
                         list_template='80px auto auto repeat(2,100px)',
                         list_hdr=['name','strain', 'intensity',  
                                   'orientation','sampling'],
                         list_fmt=['{}'] + 2*['{}'] + 2*['{}']) 
        self.pgm_options = _uiconf.get('options')   
        # output area for plots
        self._outgr = ipy.Output(layout=ipy.Layout(width='100%', border='none'))
        self.uiparam['title'] = header
        self.uiparam['list_title'] = 'Loaded data'
        self.uiparam['add_button_label'] = 'Load'
        self.uiparam['add_button_icon'] = 'file'
        self.uiparam['add_button_width'] = '80px' 
        self.uiparam['add_label'] = 'Data name'
              
        # select options
        self._options['sampling'] = _uiconf.data.listkeys('sampling')
        self._options['geometry'] = _uiconf.data.listkeys('geometry')
        # title
        self.uiparam['title'] = header

    def _init_values(self, **kwargs):
        super()._init_values()
        self._keys_conf = ['nrec']
                
    def _create_widgets(self, **kwargs):
        """Create widgets registered in self._widgets."""
        # file input
        fi = FileInput(name='strain', 
                       file=self._values['input']['strain'],
                       path=self.wk.full_path('data', as_posix=True),
                       fileonly=True,
                       label='Strain',
                       tooltip='File with strain data.',
                       width_label=80)
        self._widgets['strain'] = fi
        fi = FileInput(name='intensity', 
                       file=self._values['input']['intensity'], 
                       path=self.wk.full_path('data', as_posix=True),
                       fileonly=True,
                       label='Intensity',
                       tooltip='File with intensity data.',
                       width_label=80)
        self._widgets['intensity'] = fi

        # sampling selection
        wdg =  SelectInput(name='sampling', label='Sampling', 
                           options = self._options['sampling'],
                           value = self._values['input']['sampling'],
                           width_label=80, width_drop=100)
        self._widgets['sampling'] = wdg
        # geometry selection
        wdg =  SelectInput(name='geometry', label='Orientation', 
                           options = self._options['geometry'],
                           value = self._values['input']['geometry'],
                           width_label=80, width_drop=100)
        self._widgets['geometry'] = wdg
        
        # number of events to plot
        wdg = ArrayInput(name='nrec', label='Events', hint='',
                         value=self._values['conf']['nrec'],
                         isInt=True, step=1000, lmin=1000,lmax=100000,
                         width_label=80, width_num=100)
        self._widgets['nrec'] = wdg 

    def _create_ui(self, **kwargs):       
        btn_run = ipy.Button(description='Show',
                                layout=ipy.Layout(width='80px'))
        btn_run.on_click(self._on_replot)
        self._buttons.append(btn_run)
        # dirty trick: this is not a button, but we want it next to the button
        self._buttons.append(self._widgets['nrec'].ui())
        layout = ipy.Layout(margin='10px 5px 5px 10px')
        box = ipy.VBox([self._widgets['strain'].ui(),
                        self._widgets['intensity'].ui(),
                        self._widgets['sampling'].ui(),
                        self._widgets['geometry'].ui()],
                        layout=layout)               
        return box


    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given ID key."""
        dlist = _uiconf.data.get_item(self.name)
        if key in dlist:   
            lst = dlist[key]
            ori = lst['geometry']
            sam = lst['sampling']
            try:
                dt = lst['fdata'].scan
                if dt['eps'] is not None:
                    epsfile = dt['epsfile']
                else:
                    epsfile='none'
                if dt['int'] is not None:
                    intfile = dt['intfile']
                else:
                    intfile='none'
            except:
                epsfile = 'unknown'
                intfile = 'unknown'
        else:    
            lst = self._values['list'][key]
            ori = lst['geometry']
            sam = lst['sampling']
            epsfile = 'unknown'
            intfile = 'unknown'  
        rec = [key, epsfile, intfile, ori, sam ]
        return rec 
            
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
        
    def add_data(self, item, data=None, update_ui=False):
        """Add a dataset to the list.
        
        At the end, calls self._call_change with the list content.
        
        Parameters
        ----------
        item : str
            ID string for new item. It must be a unique ID string.
        data : dict
            New keys with values to be added to the table.
            If None, then values are taken from corresponding widgets.
        update_ui : bool
            If True, call self.update_widgets() to redraw the UI.
        """
        try:
            super().add_data(item, data=data, update_ui=False)
            _uiconf.data.load('data', item=item)
            if update_ui:
                self.update_widgets()
        except Exception as e:
            if item in self._values['list']:
                self._delete(item)
            self.error(e)
    
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
        self.message('')
        self.error('')
        self._outgr.clear_output()
        # update values from widgets
        self.update_values()
        # get command parameters
        par = _uiconf.get(self.name)
        _uiconf.data.reload_all()
        if not _uiconf.is_ready():
            self._out.clear_output()
            self.error('Input data not ready.')
            return
        # set attenuation
        att = _uiconf.data.get_attenuation()
        comm.set_attenuation(att)       
        # collect simulated pseudo-strains and experimental data to plot 
        dlist = _uiconf.data.get_item(self.name)
        expdata = {}
        simdata = {}
        for name in dlist:
            scan = _uiconf.data.get_scan(name)
            x = scan['eps'][:,0]        
            scan_range = [min(x), max(x), 2*len(x)+1]
            fname = self._get_output_filename(name=scan['epsfile'])       
            nrec = par['nrec']
            save = self.pgm_options['save']
            # set geometry and sampling for given scan
            geometry = {k:scan[k] for k in Geometry.input_keys}
            comm.set_geometry(geometry)
            comm.set_sampling(scan['sampling'])
            with self._msg:
                res = comm.report_pseudo_strains(scan_range, fname, 
                                                 nev=nrec,
                                                 intensity=True,
                                                 inline=True,
                                                 plot=False, 
                                                 save=save)
            expdata[name] = scan
            simdata[name] = res
        # do plot
        with self._outgr:
            self._plot_comparison(simdata, expdata)
      
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg, **kwargs)
        display(self._outgr)
        
        
class UI_distribution(UI_base_list):  
    """Define a free model distribution by a set of nodes and interpolation."""
    
    _dtypes = {'imodel': {'title': 'Intensity model'},
               'emodel': {'title': 'Strain model'}
               }
    
    def __init__(self, dtype='imodel'):
        if not dtype in UI_distribution._dtypes:
            raise Exception('Unkown model type: {}'.format(dtype))
        self.dtype = dtype
        super().__init__(dtype, ['dist','scale','interp'],
                         list_template='120px repeat(4,auto) ',
                         list_hdr=['name','nodes', 'range', 'variables', 
                                   'interpolation'],
                         list_fmt=['{}','{:d}', '({:g},{:g})','{:d}', '{}'] )       
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none'))         
        # title
        self.uiparam['title'] = UI_distribution._dtypes[self.dtype]['title']
        self.uiparam['list_title'] = self.uiparam['title'] + ' list'
        self.uiparam['add_button_label'] = 'Add'
        self.uiparam['add_button_icon'] = None
        self.uiparam['add_button_width'] = '50px'
        self.uiparam['add_label'] = 'Model name'
        self._was_displayed = False
        
    def _create_widgets(self, **kwargs): 

        # distribution table
        self._widgets['dist'] = DistTable(name='dist', label='',
                                          value=self._values['input']['dist'],
                                          border='1px solid')
        # scale table
        self._widgets['scale'] = ScaleTable(name='scale', label='Scaling',
                                            value=self._values['input']['scale'])
        # interpolation options
        intps = mc.intpmodels
        wdg =  SelectInput(name='interp', label='Interpolation', 
                           options = intps,
                           value = self._values['input']['interp'],
                           width_label=80, width_drop=100)        
        self._widgets['interp'] = wdg

    def _create_ui(self, **kwargs): 
        left = ipy.VBox([self._widgets['dist'].ui(), 
                         self._widgets['interp'].ui()]
                        )
        right = self._widgets['scale'].ui()
        box = ipy.HBox([left, right],
                        layout=ipy.Layout(width='100%', border='none', 
                                          padding='0px'))
        
        return box        

    def _get_display_record(self, key):
        """Return a list of values to be shown on the list for given ID key."""
        dlist = self._values['list']
        if key in dlist:   
            lst = dlist[key]
            x = lst['dist']['x']
            fitx = lst['dist']['fitx']
            fity = lst['dist']['fity']
            fs = lst['scale']['fit']
            nodes = len(x)
            rang = [x[0], x[nodes-1]]
            nfree = np.sum(fitx)
            nfree += np.sum(fity)
            nfree += np.sum(fs)
            intp = lst['interp']
        else:    
            nodes = 0
            rang = [0, 0]
            nfree = 0
        rec = [key, nodes, rang, nfree, intp]
        return rec 

    def refresh(self):
        """Redraw the model table."""
        if not self._was_displayed:
            self._widgets['dist'].redraw()
            self._widgets['scale'].redraw()
            self._was_displayed = True
    
    def show(self, err=None, msg=None, **kwargs): 
        """Display the input collection."""
        super().show(err=err, msg=msg)
        self._widgets['dist'].redraw()
        #self._widgets['dist'].sheet.layout.height='400px'
        #print(self._widgets['dist'].sheet.layout)

#%% Command collections


class UI_plot_scene(UI_base):
    """Plot scene with sample geometry and sampling events.
    
    A dialog for plotting of the sample contours with scattering geometry 
    and sampling events in several projections. 
    
    Parameters
    ----------
    rang : int
        Display range in mm. 
    proj : int
        Initially selected projection (0 ..2) 
    nrec : int
        Number of sampling events to show.
        
    Example
    -------
    
    `plot = UI_plot_scene('plot')`
    
    `plot.show()`

    """
    
    def __init__(self, name, rang=16, proj=1, nrec=3000):
        super().__init__(name, ['sampling','geometry','nrec','proj','rang'],
                         rang=rang, proj=proj, nrec=nrec)    
        self.pgm_options = _uiconf.get('options')

        # output area
        layout = ipy.Layout(width='100%', border='none', 
                            margin='0px 0px 0px 20px')
        self._out = ipy.Output(layout=layout) 
        # options linked to data lists
        self._options['sampling'] = _uiconf.data.listkeys('sampling')
        self._options['geometry'] = _uiconf.data.listkeys('geometry')
        # title
        self.uiparam['title'] = 'Plot scene'        
                
    def _create_widgets(self):         
        # sampling selection
        wdg = SelectInput(name='sampling', label='Sampling',
                          options=self._options['sampling'], 
                          value=self._values['sampling'],
                          width_label=100, width_drop=100)
        #wdg = create_select(name='sampling', label='Sampling', 
        #                    options=self._options['sampling'], 
        #                    width_label=100, width_drop=100)
        self._widgets['sampling'] = wdg
        # orientation selection
        wdg = SelectInput(name='geometry', label='Orientation',
                          options=self._options['geometry'], 
                          value=self._values['geometry'],
                          width_label=100, width_drop=100)
        #wdg =  create_select(name='geometry', label='Orientation', 
        #                     options=self._options['geometry'], 
        #                     width_label=100, width_drop=100)
        self._widgets['geometry'] = wdg
        # number of events to plot        
        wdg = ArrayInput(name='nrec', label='Events',
                         value=self._values['nrec'], hint='',
                         isInt=True, step=1000, lmin=1000,lmax=100000,
                         width_label=100, width_num=100)        
        self._widgets['nrec'] = wdg 
        # projection plane
        wdg = SelectInput(name='proj', label='Projection: ',
                          options=[('z,y',0),('x,z',1),('x,y',2)], 
                          value=self._values['proj'],
                          width_label=100, width_drop=100)
        self._widgets['proj'] = wdg 
        # plot range
        wdg = ipy.IntSlider(value=self._values['rang'], min=3, max=50, step=1, 
                            description='Range: ')
        self._widgets['rang'] = wdg 
        
    def _create_ui(self, **kwargs):        
        btn_replot = ipy.Button(description='Replot',
                                layout=ipy.Layout(width='80px'))
        btn_replot.on_click(self._on_replot)
        self._buttons.extend([btn_replot, self._widgets['rang']])
        
        box1 = ipy.VBox([self._widgets['sampling'].ui(), 
                         self._widgets['geometry'].ui(), 
                         self._widgets['proj'].ui(),
                         self._widgets['nrec'].ui()],
                        layout=ipy.Layout(width='30em'))
        layout = ipy.Layout(justify_content='center')
        box = ipy.HBox([box1, self._out],layout=layout)
        #box = ipy.VBox([self._widgets['rang'], hbox])
        return box
    
    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        pfx = self.pgm_options['prefix']
        ori = self._widgets['geometry'].value
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
        self._out.clear_output()
        # update values from widgets
        self.update_values()
        # get command parameters
        par = _uiconf.get(self.name)
        _uiconf.data.reload('sampling', item=par['sampling'])
        if not _uiconf.is_ready():
            self._out.clear_output()
            self.error('Input data not ready.')
            return      
        # set selected sampling and geometry  
        g = _uiconf.data.get_item('geometry',item=par['geometry'])
        s = _uiconf.data.get_item('sampling',item=par['sampling'])['fdata']
        comm.set_geometry(g)
        comm.set_sampling(s)
        # do plot
        save = self.pgm_options['save']
        if save:
            fname = self._get_output_filename(ext='.png')
        else:
            fname = ''
        with self._out:
            comm.plot_scene(par['nrec'], 
                            filename=fname, 
                            rang=2*[par['rang']],
                            proj=par['proj'])


class UI_resolution(UI_base):
    """Calculate resolution effects.
        
    A dialog for calculation of resolution effects: pseudo-strain, spatial
    resolution, sampling centre of mass, etc.
    
    Parameters
    ----------
    rang : int
        Display scan range in mm. 
    steps : int
        number of steps for the scan
    nrec : int
        Number of sampling events to show.
        
    Example
    -------
    
    `resol = UI_resolution('resolution')`
    
    `resol.show()`
    
    
    """
    
    def __init__(self, name, rang=[-10, 10], steps=21, nrec=3000):
        super().__init__(name, 
                         ['sampling','geometry','nrec','strain','resolution',
                          'rang','steps'],
                         rang=rang, steps=steps, nrec=nrec) 
        self.pgm_options = _uiconf.get('options')        
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none'))         
        # select options
        self._options['sampling'] = _uiconf.data.listkeys('sampling')
        self._options['geometry'] = _uiconf.data.listkeys('geometry')
        # title
        self.uiparam['title'] = 'Calculate resolution effects'         

    def _create_widgets(self, **kwargs): 

        # sampling selection
        wdg = SelectInput(name='sampling', label='Sampling',
                          options=self._options['sampling'], 
                          value=self._values['sampling'],
                          width_label=80, width_drop=100)
        #wdg = create_select(name='sampling', label='Sampling', 
        #                    options=self._options['sampling'], 
        #                    width_label=80, width_drop=100)
        self._widgets['sampling'] = wdg
        # oreiantation selection
        wdg = SelectInput(name='geometry', label='Orientation',
                          options=self._options['geometry'], 
                          value=self._values['geometry'],
                          width_label=80, width_drop=100)
        #wdg =  create_select(name='geometry', label='Orientation', 
        #                     options=self._options['geometry'], 
        #                     width_label=80, width_drop=100)
        self._widgets['geometry'] = wdg
        # number of events to plot
        wdg = ArrayInput(name='nrec', label='Events',
                         value=self._values['nrec'], hint='',
                         isInt=True, step=1000, lmin=1000,lmax=100000,
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
        
        btn_run = ipy.Button(description='Run',
                                layout=ipy.Layout(width='80px'))
        btn_run.on_click(self._on_replot)
        self._buttons.append(btn_run)
        
        layout = ipy.Layout(margin='10px 5px 5px 10px')
        box1 = ipy.VBox([self._widgets['sampling'].ui(), 
                         self._widgets['geometry'].ui(), 
                         self._widgets['nrec'].ui()],
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
        ori = self._widgets['geometry'].value
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
        self._out.clear_output()
        # update values from widgets
        self.update_values()
        
        # get command parameters
        par = _uiconf.get(self.name)
        _uiconf.data.reload('sampling', item=par['sampling'])
        _uiconf.data.reload('attenuation')
        if not _uiconf.is_ready():
            self._out.clear_output()
            self.error('Input data not ready.')
            return

        # set selected sampling and geometry
        g = _uiconf.data.get_item('geometry',item=par['geometry'])
        s = _uiconf.data.get_item('sampling',item=par['sampling'])['fdata']
        comm.set_geometry(g)
        comm.set_sampling(s)
        # set attenuation
        att = _uiconf.data.get_attenuation()
        comm.set_attenuation(att)

        # do plot
        nrec = par['nrec']
        fname = self._get_output_filename()
        rang = list(par['rang'])
        nstp = par['steps']
        scan_range = rang + [nstp]
        save = self.pgm_options['save']
        with self._out:
            if par['strain']:
                comm.report_pseudo_strains(scan_range, fname, 
                                           nev=nrec,
                                           intensity=True,
                                           inline=True,
                                           plot=True, 
                                           save=save)
            if par['resolution']:    
                comm.report_resolution(scan_range, fname, 
                                       nev=nrec,
                                       cog=True,
                                       inline=True,
                                       plot=True, 
                                       save=save)
    


class UI_fit_imodel(UI_base):
    """Fit intensity model.
    
    Parameters
    ----------
    rang : int
        Display scan range in mm. 
    steps : int
        number of steps for the scan
    nrec : int
        Number of sampling events to show.
        
    Example
    -------
    
    `resol = UI_resolution('resolution')`
    
    `resol.show()`
    
    
    """
    
    def __init__(self, name):
        super().__init__(name, 
                         ['data','model','nrec','npts', 'fit']) 
        self.pgm_options = _uiconf.get('options')        
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none'))         
        # select options
        self._options['data'] = _uiconf.data.listkeys('data')
        self._options['model'] = _uiconf.data.listkeys('imodel')
        # title
        self.uiparam['title'] = 'Fit intensity model'

    def _create_widgets(self, **kwargs): 

        # sampling selection
        wdg = SelectInput(name='data', label='Data',
                          options=self._options['data'], 
                          value=self._values['data'],
                          width_label=80, width_drop=100)
        self._widgets['data'] = wdg
        # oreiantation selection
        wdg = SelectInput(name='model', label='Model',
                          options=self._options['model'], 
                          value=self._values['model'],
                          width_label=80, width_drop=100)
        #wdg =  create_select(name='geometry', label='Orientation', 
        #                     options=self._options['geometry'], 
        #                     width_label=80, width_drop=100)
        self._widgets['model'] = wdg
        # number of events to plot
        wdg = ArrayInput(name='nrec', label='Events',
                         value=self._values['nrec'], 
                         hint='',
                         isInt=True, step=1000, lmin=1000,lmax=100000,
                         width_label=80, width_num=100)   
        self._widgets['nrec'] = wdg 
        # number of points on model curve
        wdg = ArrayInput(name='npts', label='Points',
                         value=self._values['npts'], 
                         hint='',
                         isInt=True, step=10, lmin=10,lmax=500,
                         width_label=80, width_num=100)   
        self._widgets['npts'] = wdg
        
        # fit control parameters
        wdg = FitControl(name='fit', value=self._values['fit'])
        self._widgets['fit'] = wdg

    def _create_ui(self, **kwargs):
        
        btn_run = ipy.Button(description='Plot',
                                layout=ipy.Layout(width='80px'))
        btn_run.on_click(self._on_replot)
        self._buttons.append(btn_run)
        
        btn_guess = ipy.Button(description='Guess',
                                layout=ipy.Layout(width='80px'))
        btn_guess.on_click(self._on_guess)
        #self._buttons.append(btn_guess)
        
        btn_fit = ipy.Button(description='Fit',
                                layout=ipy.Layout(width='80px'))
        btn_fit.on_click(self._on_fit)
        self._buttons.append(btn_fit)
        
        layout = ipy.Layout(margin='10px 10px 5px 10px')
        box1 = ipy.VBox([self._widgets['data'].ui(), 
                         self._widgets['model'].ui(), 
                         self._widgets['nrec'].ui(),
                         self._widgets['npts'].ui()],
                         layout=layout)               
        box2 = ipy.VBox([self._widgets['fit'].ui()],
                         layout=layout) 
        #box2 = ipy.VBox([self._widgets['strain'], 
        #                 self._widgets['resolution']],
        #                 layout=layout)
        #layout=ipy.Layout(margin='10px 5px 5px 10px', min_width='30%')
        #box3 = ipy.VBox([self._widgets['rang'].ui(), 
        #                 self._widgets['steps'].ui()],
        #                 layout=layout)
        
        #style = {'margin': '5px 5px 5px 50px'}
        hbox = ipy.HBox([box1,box2])
        box = ipy.VBox([hbox, self._out])
        return box
            
    def _get_output_filename(self, ext=''):
        """Generate base output filename."""
        pfx = self.pgm_options['prefix']
        dname = self._widgets['data'].value
        
        if pfx:
            base = '{}_{}_imodel'.format(pfx, dname)
        else:
            base = '{}_imodel'.format(dname)
        if ext:
            fname = base + ext
        else:
            fname = base
        return fname
    
    def _prepare(self):
        """Make all settings necessary to run commands."""
        self.message('')
        self.error('')
        self._out.clear_output()
        # update values from widgets
        self.update_values()
        # reload all input data if needed
        _uiconf.data.reload_all()
        # get command parameters
        par = _uiconf.get(self.name)
        # selected scan data 
        scan = _uiconf.data.get_scan(par['data'])
        # selected model
        model = _uiconf.data.get_item('imodel',item=par['model'])
        # scale parameters
        scale = model['scale']
        
        
        # set scan parameters and sampling
        comm.set_scan(scan)        
        # set attenuation
        att = _uiconf.data.get_attenuation()
        comm.set_attenuation(att)        
        # create ifit object
        ifit = comm.define_ifit(scan, model['dist'], par['nrec'], 
                                ndim=par['npts'])        
        # define scaling parameters and interpolation model
        ifit.defScaling(scale['values'], scale['fit'], 
                        minval=[0., 0., -np.inf])
        ifit.setInterpModel(model['interp'])
        return ifit
    
    def _on_replot(self,b):
        ifit = self._prepare()        
        with self._out:
            # create smeared curve: run without fitting, maxiter=0
            comm.run_fit(ifit, maxiter=0)
            # plot results
            comm.report_fit(ifit, '')
        
    def _on_guess(self,b):
        ifit = self._prepare()    
        par = _uiconf.get(self.name)
        with self._out:
            # create smeared curve: run without fitting, maxiter=0
            comm.run_fit(ifit, maxiter=par['fit']['maxiter'], guess=True)
            # plot results
            comm.report_fit(ifit, '')
            
    def _on_fit(self,b):
        ifit = self._prepare()  
        par = _uiconf.get(self.name)
        with self._out:
            # create smeared curve: run without fitting, maxiter=0
            comm.run_fit(ifit, 
                         maxiter=par['fit']['maxiter'], 
                         loops=par['fit']['loops'],
                         areg=par['fit']['areg'],
                         guess=par['fit']['guess'],
                         bootstrap=par['fit']['loops']>2)
            # plot results
            save = self.pgm_options['save']
            if save:
                fname = self._get_output_filename()
            else:
                fname = ''
            comm.report_fit(ifit, fname)
            
            
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
                      'data',
                      'imodel',
                      'fit_imodel']
    def __init__(self):    
        self._last_input_file = 'input.json'
        self.ui = {}
        self.wk = dataio.workspace()
        self._msg = ipy.Output(layout=ipy.Layout(width='100%'))
        self._err = ipy.Output(layout=ipy.Layout(width='100%'))
        # setup dialog
        self._add_input_ui(UI_workspace('workspace'))
        self._add_input_ui(UI_options('options'))
        self._add_input_ui(UI_shape('shape'))
        self._add_input_ui(UI_geometry('geometry'))
        self._add_input_ui(UI_sampling('sampling'))
        self._add_input_ui(UI_attenuation('attenuation'))
        
        # execution dialogs
        self._add_input_ui(UI_plot_scene('scene'))
        self._add_input_ui(UI_resolution('resolution')) 
        self._add_input_ui(UI_data('data'))
        self._add_input_ui(UI_distribution('imodel'))
        self._add_input_ui(UI_fit_imodel('fit_imodel'))
        
        # set event handlers
        self.ui['shape'].set_on_change(self._change)
        self.ui['geometry'].set_on_change(self._change)
        self.ui['sampling'].set_on_change(self._change)
        self.ui['data'].set_on_change(self._change)
        self.ui['imodel'].set_on_change(self._change)

    
    def _add_input_ui(self, ui):
        if ui.name in UI._registered_ui:
            self.ui[ui.name] = ui
        else:
            msg = 'Attempt to add non-registered ui: {}'
            raise Exception(msg.format(ui.name))
        
    def _tab_observe(self, widget):
        idx = widget['new']
        if idx==6:
            self.ui['imodel'].refresh()
        
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
        tabs_data['imodel'] = {'title':'Intensity',
                                  'ui':[self.ui['imodel'],
                                        self.ui['fit_imodel']]}
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
        
        tab.observe(self._tab_observe, names='selected_index')
    
    def _on_save_input(self,b):
        s = choose_file_save(initialdir=self.wk.path('work').as_posix(), 
                        initialfile=self._last_input_file,
                        filetypes=(('Setup files','*.inp')))
        if s:
            self.save(filename=s)
            p = _Path(s)
            self._last_input_file = p.name
            
    def _on_load_input(self,b):
        s = choose_file(initialdir=self.wk.path('work').as_posix(), 
                        initialfile=self._last_input_file,
                        filetypes=(('Setup files','*.inp')))
        if s:
            self.load(filename=s)
            p = _Path(s)
            self._last_input_file = p.name
            
    def _change(self, obj, **kwargs):
        if obj.name == 'shape':
            if 'shape' in kwargs:
                if self.ui['shape'].select != 'File':
                    comp = self.ui['shape'].create_shape()
                    comm.set_shape(comp)
            else:
                if self.ui['shape'].select != 'File':
                    comm.set_shape(None,**kwargs)
        # update selection lists: geometry and sampling
        self._update_options(obj.name)

    def _update_options(self, name):
        """Update selection lists: geometry and sampling."""
        if name in ['geometry', 'sampling']:
            for key in ['scene','resolution','data']:
                self.ui[key].update_options(name)
        # update other selection lists
        elif name == 'imodel':
            self.ui['fit_imodel'].update_options('model')
        elif name == 'data':
            self.ui['fit_imodel'].update_options('data')

    def save(self, filename=''):
        """Save input data in JSON format.""" 
        try:
            for key in self.ui: 
                self.ui[key].update_values()
            out = _uiconf.export_json()
            if filename:
                with open(filename, 'w') as f:
                    f.write(out)
                    f.close()
            else:
                with self._msg:
                    print(out)
        except Exception as e:
            with self._err:
                print(e)
            raise e
        
    def load(self, filename=''):
        """Load input data in JSON format."""
        try:
            # load dict from file
            with open(filename, 'r') as f:
                lines=f.readlines()
                f.close()
            # import into global input data
            _uiconf.import_json(lines, reload=True)
            if not _uiconf.is_ready():
                print('Loading of program settings failed.')
                return
            # Notification that data lists may have changed 
            # -> update dependent options
            for key in ['geometry', 'sampling', 'data', 'imodel']:
                self._update_options(key)
            # re-assign global input to the ui components
            for key in self.ui:
                try:
                    self.ui[key].assign()
                    self.ui[key].update_widgets() 
                    # self.ui[key].notify()
                except Exception as e:
                    with self._err:
                        print('Problem with settig value {}:'.format(key))
                        print(e)
                    raise e
        except Exception as e:
            with self._err:
                print(e)
            raise e

