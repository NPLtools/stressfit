# -*- coding: utf-8 -*-
"""
Classes and functions for creating user interface with ipywidgets.

Employs tkinter for file open dialogs.

Created on Tue Oct  5 11:38:15 2021
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import abc
import numpy as np
import ipywidgets as ipy
import json
from traitlets import traitlets
from pathlib import Path as _Path
from IPython.display import display
import stressfit.commands as comm
import stressfit.shapes as shapes
import stressfit.dataio as dataio
deg = np.pi/180.

def _has_keys(args, keys):
    return all (k in args for k in keys)

def choose_path(initialdir=None):
    from tkinter import Tk, filedialog
    root = Tk()
    root.withdraw() # Hides small tkinter window.
    root.attributes('-topmost', True)
    open_file = filedialog.askdirectory(initialdir=initialdir)
    return open_file

def choose_file(initialdir=None, initialfile=None):
    from tkinter import Tk, filedialog
    root = Tk()
    root.withdraw() # Hides small tkinter window.
    root.attributes('-topmost', True)
    open_file = filedialog.askopenfilename(initialdir=initialdir, 
                                           initialfile=initialfile)
    return open_file

#%% Input components

def create_input_float(label='Float number',value=0.0, 
                       lmin=None, lmax=None, step=0.1,
                       width_label=100, width_num=100, unit='px'):
    width = '{}{}'.format(width_label+width_num, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = ipy.BoundedFloatText(description=label, value=value,
                                min=lmin, max=lmax, step=step, 
                                style=label_style,
                                layout=ipy.Layout(width=width))
    return inp

def create_input_int(label='Integer',value=3000, 
                       lmin=None, lmax=None, step=1,
                       width_label=100, width_num=100, unit='px'):
    width = '{}{}'.format(width_label+width_num, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = ipy.BoundedIntText(description=label, value=value,
                                min=lmin, max=lmax, step=step, 
                                style=label_style,
                                layout=ipy.Layout(width=width))
    return inp


def create_header(text, color='#35446B', size='+1'):
    """Return header as HTML widget."""
    font_fmt = "<b><font color='{}' size='{}'>{}</font></b>"
    html = font_fmt.format(color,size,text)
    lbl = ipy.HTML(value=html)
    return lbl

   
def create_select(options, label, value=None, index=0, width_label=100, 
                  width_drop=60, unit='px'):
    """
    Create Dropdown widget for given options.

    Parameters
    ----------
    options : list of (name, value)
        Selection options.
    label : str
        Label text on the left.
    value : 
        Initial selected value. Overrides index.
    index : int
        Initial selection index.
    width_label : int
        Width of label text in given units.
    width_drop :int
        Width of dropbox widget in given units.
    unit : str
        Width units.

    Returns
    -------
    sel : widgets.Dropdown

    """
    if value is not None:
        vini=value
    elif len(options)>index:
        vini = options[index][1]
    else:
        vini=None
    max_width = '{}{}'.format(width_label+width_drop, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    sel = ipy.Dropdown(description=label,options=options,value=vini,
                           layout=ipy.Layout(max_width=max_width),
                           style=label_style)
    return sel


class ValueButton(ipy.Button):
    """Button with value trait."""
    def __init__(self, value='', *args, **kwargs):
        super(ValueButton, self).__init__(*args, **kwargs)
        self.add_traits(value=traitlets.Any(value))

class DirInput():
    """Directory input as Text widget, open button and tooltip.
    
    Parameters
    ----------
    path : str
        Path string
    label : str
        Label to appear before the text input.
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    output : `ipywidgets.Output`
        Optional output widget is used to show error messages.
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    
    """
    def __init__(self, path, label='', tooltip='', output=None,
                 width_label=150, width_button=50):
        self.path = path
        self.label = label
        self.tooltip = tooltip
        self.output = output
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = ipy.Text(self.path,layout=ipy.Layout(width='100%'))
        self.btn = ValueButton(description='', value=self.txt.value, 
                               layout=ipy.Layout(width=w_btn),
                               tooltip=self.tooltip, icon='folder-open')
        # link button and text
        self.dl = ipy.link((self.txt, 'value'), (self.btn, 'value'))
        # add events
        self.btn.on_click(self._on_button)
        self.txt.observe(self._on_text,'value', type='change')

    def _on_button(self,ex):
        """Open path dialog and save result."""
        s = choose_path(initialdir=ex.value)
        if s:
            ex.value=s
            self.path = s
            if self.output:
                self.output.clear_output()
    
    def _on_text(self, change):
        """Text change - save result."""
        s = change['new']
        if s:
            self.path = s
            if self.output:
                self.output.clear_output()
    
    def get_value(self):
        """Return path as pathlib.Path object."""
        return _Path(self.path)
    
    def set_value(self, value):
        if isinstance(value, _Path):
            self.txt.value = value.as_posix()
        else:
            self.txt.value = str(value)
    
    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui
   

class FileInput():
    """File input as Text widget, open button and tooltip.
    
    Parameters
    ----------
    file : str
        Filename - full path or basename.
    path : str
        Optional directory name. 
    label : str
        Label to appear before the text input.
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    output : `ipywidgets.Output`
        Optional output widget is used to show error messages.
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    
    """
    def __init__(self, file, path=None, label='', tooltip='', output=None,
                 width_label=150, width_button=50):
        
        if path is None:
            self.path = dataio.get_paths()['data']
        else:
            self.path = path
        self._set_file(file)
        self.label = label
        self.tooltip = tooltip
        self.output = output
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = ipy.Text(self.file,layout=ipy.Layout(width='100%'))
        self.btn = ValueButton(description='', value=self.txt.value, 
                               layout=ipy.Layout(width=w_btn),
                               tooltip=self.tooltip, icon='folder-open')
        # add events
        self.btn.on_click(self._on_button)

    def _set_file(self, file):
        """Set new file name.
        
        If file is relative, join it with path info. 
        
        If file is absolute, update the path info.                
        """
        f = _Path(file)
        if f.is_absolute():
            p = f.parent
            self.fullname = f
            self.file = f.name
            self.path = p.as_posix()
        else:
            p = _Path(self.path)
            self.fullname = p.joinpath(file)
            try:
                self.file = self.fullname.relative_to(p).as_posix()
            except:
                self.file = file # this should not happen                   
        
    def _on_button(self,ex):
        """Open path dialog and save result."""
        s = choose_file(initialdir=self.fullname.parent.as_posix(), 
                        initialfile=self.fullname.name)
        if s:
            self._set_file(s)
            ex.value=s
            self.txt.value = self.file
            if self.output:
                self.output.clear_output()
                
    def get_value(self):
        """Return file name as dict with file and path items."""
        return {'path': self.path, 'file': self.file}    
    
    def set_value(self, value):
        """Set file name from dict with file and path items.
        
        The path item is optional. 
        """
        if 'path' in value and value['path'] is not None:
            self.path = value['path']
        self._set_file(value['file'])

    
    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui
     

class PathInput():
    """Path input as Text widget, open button and tooltip.
    
    Saves result in the dict "dataset" under given key passed as arguments.
    
    Parameters
    ----------
    dataset : dict
        Dictionary where to store result.
    key: str
        Key under which the result is stored.
    label : str
        Label to appear before the text input.
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    file : bool
        If true, input should be a file name. 
    output : `ipywidgets.Output`
        Optional output widget is used to show error messages.
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    callback : func(str)
        call-back function invoked after by button click.
    kwargs : dict
        other arguments passed to Text constructor.
    
    """
    def __init__(self, dataset, key, label='', tooltip='', output=None,  
                 width_label=150, width_button=50, callback=None, **kwargs):
        self.dataset = dataset
        self.key = key
        self.label = label
        self.tooltip = tooltip
        self.output = output
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = ipy.Text(self.dataset[self.key], 
                            layout=ipy.Layout(width='100%'), **kwargs)
        self.btn = ValueButton(description='', value=self.txt.value, 
                               layout=ipy.Layout(width=w_btn),
                               tooltip=self.tooltip, icon='folder-open')
        # link button and text
        self.dl = ipy.link((self.txt, 'value'), (self.btn, 'value'))
        # add events
        self.btn.on_click(self.on_button)
        self.txt.observe(self.on_text,'value', type='change')
        self.callback = callback

    def on_button(self,ex):
        """Open file or path dialog and save result."""
        s = choose_path(initialdir=ex.value)
        if s:
            ex.value=s
            self.dataset[self.key] = s
            if self.output:
                self.output.clear_output()
            if self.callback:
                self.callback(s)
                
    
    def on_text(self, change):
        """Text change - save result."""
        s = change['new']
        if s:
            self.dataset[self.key] = s
            if self.output:
                self.output.clear_output()
    
    def set_value(self, value:str):
        """Set new path name and update widgets.
        """
        self.txt.value = value
        self.btn.value = value
      
    def get_value(self):
        return self.txt.value
    
    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui
    
class ArrayInput():
    """Array input set.
    
    Input set for a numerical array formatted as GridBox. Includes
    label, a row of numeric inputs and tooltip text.
    Scalar value is also accepted.
    
    Parameters
    ----------
    value: array_like
        Input array. The legth defines number of input widgets.
    label: str
        Label text on the left.
    hint: str
        Tooltip text on the right.
    step: float
        For float input,  defines step. 
    isInt: bool
        If true, the aray eleents are treated as integers.
    width_label : int
        Width of label text in given units.
    width_num :int
        Width of numeric input widget in given units.
    unit : str
        Width units.
    """
    def __init__(self, value=[0.0,0.0,0.0], label='array', 
                 hint='define an array', step=0.1, isInt=False,
                 width_label=150, width_num=80, unit='px'):
        if isinstance(value, (float,int)):
            self.value = [value]
            self._is_scalar = True
        else:
            self.value = value
            self._is_scalar = False
        self.label = label
        self.hint = hint
        self.isInt = isInt
        self.inputs = []
        x = [width_label,unit,len(self.value), width_num, unit]
        self._tpl = '{}{} repeat({},{}{}) auto'.format(*x)
        for j in range(len(self.value)):
            w = '{}{}'.format(width_num,unit)
            layout=ipy.Layout(width=w)
            if self.isInt:
                b = ipy.IntText(description='',value=int(self.value[j]),
                                    step=1, layout=layout)
            else:
                b = ipy.FloatText(description='',value=float(self.value[j]),
                                      step=step, layout=layout)
            self.inputs.append(b)
    
    def get_value(self):
        """Get values from the widgets."""
        if self._is_scalar:
            arr = self.inputs[0].value
        else:
            arr = np.zeros(len(self.inputs))
            for i in range(len(self.inputs)):
                arr[i] = self.inputs[i].value
        return arr
    
    def set_value(self, data):
        """Set values to the widgets."""
        if self._is_scalar:
            if not (isinstance(data, float) or isinstance(data, int)):
                raise Exception('Required scalar value, not {}'.format(data))
            self.inputs[0].value = data
        else:
            nd = len(data)
            ni = len(self.inputs)
            if len(data) != len(self.inputs):
                raise Exception('Array lengths do not match: {},{}'.format(nd,ni))
            for i in range(ni):
                self.inputs[i].value = data[i]
    
    def ui(self):
        wdg = []
        grid_layout = ipy.Layout(grid_template_columns=self._tpl, 
                                     grid_gap='5px 5px')
        lbl = ipy.HTML('{}'.format(self.label))
        hint = ipy.Label(self.hint)
        wdg = []
        wdg.append(lbl)
        wdg.extend(self.inputs)
        wdg.append(hint)
        box = ipy.GridBox(wdg, layout=grid_layout)          
        return box   


class SelectInput():
    """Dropdown input set.
    
    Input set for an option list formatted as GridBox. Includes
    label, a Dropdown and tooltip text.
    
    Parameters
    ----------
    options : list of (name, value)
        Selection options.
    label : str
        Label text on the left.
    hint: str
        Tooltip text on the right.
    value : 
        Initial selected value. Overrides index.
    index : int
        Initial selection index.
    width_label : int
        Width of label text in given units.
    width_drop :int
        Width of dropbox widget in given units.
    unit : str
        Width units.
    """
    def __init__(self, options, label, hint='', value=None, index=0, 
                 width_label=150, width_drop=60, unit='px'):
        self.options = options
        self.label = label
        self.hint = hint
        self.value = value
        self.index = index
        x = [width_label,unit, width_drop, unit]
        self._tpl = '{}{} {}{} auto'.format(*x)
        vini = None
        if value is not None:
            vini=value
        elif index>-1 and len(options)>index:
            if isinstance(options[index], (float, int, str)):
                # options is list of strings
                vini = options[index]
            elif len(options[index])>1:
                # options is list of (name, value)
                vini = options[index][1]

        w = '{}{}'.format(width_drop,unit)
        layout = ipy.Layout(width=w)
        self.input = ipy.Dropdown(description='', options=options, 
                                      value=vini, layout=layout)
    
    def get_value(self):
        """Get selected value."""
        return self.input.value

    def set_value(self, data):
        """Set value to the widget."""
        try:
            self.input.value = data
        except Exception as e:
            msg = 'Cannot set value to options: {}\n{}'
            raise Exception(msg.format(data,e))

    def set_index(self, index):
        """Set value by index to the widget."""
        if index >= len(self.options):
            msg = 'Index exceeds options length: {}>{}'
            raise Exception(msg.format(index,len(self.options)))
        self.input.index = index
       
    def get_index(self):
        """Get selected index."""
        return self.input.index
    
    def ui(self):
        wdg = []
        grid_layout = ipy.Layout(grid_template_columns=self._tpl, 
                                     grid_gap='5px 5px')
        lbl = ipy.HTML('{}'.format(self.label))
        hint = ipy.Label(self.hint)
        wdg = []
        wdg.append(lbl)
        wdg.append(self.input)
        wdg.append(hint)
        box = ipy.GridBox(wdg, layout=grid_layout)          
        return box  

#%% Input collections



class UI_base:
    """Base abstract class for input collection."""
    
    def _err_ni(src):
        msg = "{}(): Subclass must implement abstract method".format(src)
        raise NotImplementedError(msg)
    
    _NON_IMPLEMENTED = "Subclass must implement abstract method"
    def __init__(self, name, keys):
        self._name = name
        self._keys = keys # list of allowed value keys
        self._widgets = {} # dict of input widgets
        self._values = {} # dict of widget values
        self._on_change = None
    
    def _check_keys(self,data):
        """Verify that data contains all required keys."""
        if not _has_keys(data,self._keys):
            msg = 'Required key(s) missing in the input data: {}'
            raise Exception(msg.format(self._keys))
            
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

    @abc.abstractmethod
    def show(self):
        """Create widgets if not done in constructor and display."""        
        UI_base._err_ni('show')

    ### Other methods
    def get_values(self):
        """Return actual input values as dict."""
        self.update_values()
        return self._values

    def set_on_change(self, on_change):
        """Set callback to be executed when the UI changes state."""
        self._on_change = on_change
        
    def on_change(self, **kwargs):
        """Call when the UI changes state."""
        if self._on_change:
            self._on_change(self, **kwargs)

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
    style_lbl = {'description_width': '150px'}
    def __init__(self, select=shapes.Plate):
        super().__init__('shape', ['shape', 'param'])   
        self.shapes = dataio.load_config('shapes')
        self.select = select
        self.options = []
        for key in self.shapes:
            self.options.append((self.shapes[key]['name'],key))
        
        self._sel = create_select(self.options, label='', width_label=0, 
                            width_drop=150, unit='px')
        self._sel.value = self.select
        self._sel.observe(self._on_change_shape, type='change')
        self._out = ipy.Output(layout=ipy.Layout(border='1px solid'))

    def _reset_widgets(self):
        """Re-fill self._widgets with input widgets for selected shape."""
        param = self.shapes[self.select]['param']
        self._widgets.clear()
        self._widgets['select']=self.select
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
                wgt = ArrayInput(value=value, label=desc, hint=hint)
            elif isinstance(value, list):
                if isinstance(value[0],float):
                    wgt = ArrayInput(value=value, label=desc, hint=hint)
                else:
                    wgt = SelectInput(options=value, label=desc, hint=hint,
                                      width_drop=80)
            else:
                raise Exception('Unknown parameter type: {}'.format(value))
            items[key] = wgt
        
        self._widgets['param'] = items

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
            self.on_change()
 
    def update_widgets(self, data):
        """Update parameter input widgets from data.

        Parameters
        ----------
        data: dict
            Input names and values 
        
        """
        self._check_keys(data)
        if data['select'] != self.select:
            self._on_change_shape({'name':'value', 'new':data['select']})
        param = data['param']
        items = self._widgets['param']
        for key in param:
            item = items[key]
            if isinstance(item, SelectInput):
                items[key].set_index(param[key].index)
            else:
                items[key].set_value(param[key].value)
        self.update_values()            
        
    def update_values(self):
        """Update input data from widgets.
        
        Input data are stored in self._values as dict.
        
        """
        self._values.clear()
        self._values['select']=self.select
        param = {}
        items = self._widgets['param']
        for key in items:
            item = items[key]
            if isinstance(item, SelectInput):
                value = item.get_index()
            else:
                value = item.get_value()
            param[key] = value
        self._values['param'] = param
    
    def show(self):   
        lbl = create_header('Sample shape')
        display(lbl)
        display(self._sel)
        display(self._out)
        self._on_change_shape({'name':'value', 'new':self.select})
    
    def create_shape(self):
        """Create selected shape instance from stressfit.shapes."""
        self.update_values()
        comp = shapes.create(self._values['select'],**self._values['param'])
        return comp

        
class UI_orientation(UI_base):
    """Dialog for sample orientation.
    
    This UI allows to set the sample and scan orientation and mainain
    a named list of such oreiantations.
    
    Each orientation is defined by four vectors:    
    
    angles : array_like(3)
        Euler YZY angles of sample orientation
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
        super().__init__('orientation', 
                         ['angles', 'scandir', 'rotctr', 'scanorig'])   
        # widgets for current dataset
        self._widgets['angles'] = ArrayInput(value=[135,0,0],label='angles',
                            hint='Sample rotation: YZY Euler angles [deg]')
        self._widgets['scandir'] = ArrayInput(value=[0,0,1],label='scandir', 
                            hint='Scan direction vector')
        self._widgets['rotctr'] = ArrayInput(value=[0,0,0],label='rotctr', 
                            hint='Rotation centre [mm]')
        self._widgets['scanorig'] = ArrayInput(value=[0,0,0],label='scanorig', 
                            hint='Scan origin [mm]')
        self._values['input'] = {}
        self._values['list'] = {}
        self._out = ipy.Output(layout=ipy.Layout(border='1px solid'))
        self._name_inp = ipy.Text(description='',value='', 
                                 layout=ipy.Layout(width='100px'))
    
    
    def _get_ori(self):
        """Retrieve orientation data from widgets."""
        vals = {}
        for key in self._widgets:
            vals[key] = self._widgets[key].get_value()
        return vals
    
    def _on_add(self,b):
        txt = self._name_inp.value
#print(txt)
        if not txt:
            print('Define a unique name')
        elif not txt or txt in self._values['list']:
            print('The orientation name "{}" is already defined'.format(txt))
        else:
            self.add_orientation(txt, update_table=True)

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
            b = ValueButton(description='',icon='trash',value=key,
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
                Values for the orientation input widgets (four vectors).
            list : dict, optional
                Orientation datasets to be added in the table.
                Each dataset contains the four vectors defining an orientation.
        
        """
        if 'input' in data:
            dset = data['input']
            self._check_keys(dset)
            for key in self._keys:
                self._widgets[key].set_value(dset[key])
        
        if 'list' in data:
            lst = data['list']
            self._values['list'].clear()
            for key in lst:
                self._check_keys(lst[key])
                self.add_orientation(lst[key])
            self._update_table()
    
    def update_values(self):
        """Update input data from widgets.
        
        Does not change the list of defined orientations.

        """ 
        self._values['input'].clear()
        self._values['input'].update(self._get_ori())

    def show(self):
        lbl = create_header('Sample orientation')
        display(lbl)
        lst = []
        vals = list(self._widgets.values())
        for i in range(len(vals)):
            lst.append(vals[i].ui())
        box = ipy.VBox(lst)
        display(box)
        btn_add = ipy.Button(description='Add', 
                             layout=ipy.Layout(width='50px'))
        btn_add.on_click(self._on_add)
        hint = ipy.Label('Enter unique name (e.g. "radial", "axial", ...)')
        display(ipy.HBox([hint, self._name_inp, btn_add]))
        display(self._out)

    def add_orientation(self, name, update_table=False):
        """Add orientation using current input values."""
        if not name in self._values['list']:
            vals = {}
            for key in self._widgets:
                vals[key] = self._widgets[key].get_value()
            self._values['list'][name] =vals
            if update_table:
                self._update_table()
            self.on_change()
        
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
            or :meth:`sampling_list` for the named list of all 
            loaded data sets.
    
    Parameters
    ----------
    file : str
        Initial sampling file name
    path : str
        Initial directory name with sampling files
    nev : int
        Initial number of sampling points to be loaded.
    load_as : str or None
        If defined and not empty, load sampling data with initial setting.
        The name of the data set is the value of the load_as parameter.
        The attempt to load the data is made at the end of :meth:`show`.
    """
    
    style_lbl = {'description_width': '100px'}
    def __init__(self, file='', path='', nev=3000, load_as=None):
        super().__init__('sampling', ['file', 'path', 'nev'])        
        self._load_as = load_as
        self._sampling = {}
        self._values['input'] = {}
        self._values['list'] = {}
        self._values['input']['file'] = file
        self._values['input']['path'] = path
        self._values['input']['nev'] = nev
        self._out = ipy.Output(layout=ipy.Layout(width='100%', 
                                                 border='1px solid')) 
        self._err = ipy.Output(layout=ipy.Layout(width='100%', 
                                                border='none')) 
        # file input
        fi = FileInput(self._values['input']['file'], 
                       path=self._values['input']['path'], 
                       label='File name',
                       tooltip='File with simulated sampling points.')
        self._widgets['file'] = fi
        # num. of events
        inp = create_input_int(label='Maximum events',
                               value=self._values['input']['nev'],
                               lmin=1000, lmax=100000, step=1000, 
                               width_label=150, width_num=100)
        self._widgets['nev'] = inp
        # ID string
        ninp = ipy.Text(description='Name',value='', 
                                 layout=ipy.Layout(width='150px'),
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
            b = ValueButton(description='', icon='trash', value=key,
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
            self._widgets['file'].set_value(data['input'])
            self._widgets['nev'].value = data['input']['nev']
            self.update_values()
     
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
        file = self._widgets['file'].get_value()
        self._values['input']['file'] = file['file']
        self._values['input']['path'] = file['path']
        self._values['input']['nev'] = self._widgets['nev'].value

    def show(self):
        lbl = create_header('Sampling points')      
       
        # load button
        btn = ValueButton(description='Load',  icon='file', 
                          layout=ipy.Layout(width='10%'))     
        btn.on_click(self._on_load_click)
        
        # collect boxes and display
        box_add = ipy.HBox([self._widgets['nev'], self._widgets['name'], btn ])
        box = ipy.VBox([lbl, self._widgets['file'].ui(), box_add])
        display(box)
        display(self._err)
        display(self._out)
        if self._load_as:
            self.add_sampling(self._load_as, self.get_values()['input'], 
                              update_table=True )
            
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
            with self._err:
                print('Provide a unique name for new sampling.')
            return
        
        self._check_keys(data)
        sampling = comm.load_sampling(**data)
        if sampling is not None:
            data['file'] = sampling.src['file']
            data['path'] = sampling.src['path']
            data['nev'] = sampling.src['nrec']
            self._widgets['file'].set_value(data)
            self._values['input'].update(data)
            self._values['list'][name]=data.copy()
            self._sampling[name]=sampling
            if update_table:
                self._update_table()   
            self.on_change()

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
        

class UI_plot_scene(UI_base):
    """
    A dialog for plotting of the sample contours with scattering geometry 
    and sampling events in several projections. 
    
    Parameters
    ----------
    sampling_list: dict
        List of loaded sampling data sets.
        as defined by the class :class:`UI_sampling` (use :meth:`get_values`)
    orientations: dict
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
    Provided that the objects `ui_s` and `ui_o` are the dialogs for sampling
    and orientation, respectively, than create and display this dialog
    as 
    
    `plot = UI_plot_scene(ui_s.sampling_list(), ui_o.get_values)`
    
    `plot.show()`
    
    
    """
    def __init__(self, sampling_list, orientations, rang=16, proj=1, nev=3000):
        super().__init__('scene', ['sampling','ori','nrec','proj','range'])    
        self.header = create_header('Plot scene')
        self.sampling_list = sampling_list
        self.orientations = orientations
        self.options_sampling = []
        self.options_ori = []
        self._update_options()
        # sampling selection
        wdg = create_select(label='Select sampling', 
                            options=self.options_sampling, 
                            width_label=100, width_drop=80)
        self._widgets['sampling'] = wdg
        # oreiantation selection
        wdg =  create_select(label='Select orientation', 
                             options=self.options_ori, 
                             width_label=130, width_drop=100)
        self._widgets['ori'] = wdg
        # number of events to plot
        wdg = create_input_int(label='Events',
                               value=nev,
                               lmin=1000, lmax=100000, step=1000, 
                               width_label=150, width_num=100)
        self._widgets['nrec'] = wdg 
        # projection plane
        wdg = create_select(label='Projection: ', 
                            options=[('z,y',0),('x,z',1),('x,y',2)],
                            value = proj,
                            width_label=100, width_drop=80)
        self._widgets['proj'] = wdg 
        # plot range
        wdg = ipy.IntSlider(value=rang, min=3, max=100, step=1, 
                            description='Range: ')
        self._widgets['range'] = wdg 
        # output area
        self._out = ipy.Output(layout=ipy.Layout(width='100%', border='none')) 
        
    def _update_options(self):
        """Update options for drop-down lists"""
        self.options_sampling.clear()
        #for i in range(len(self.sampling_list)):
        for key in self.sampling_list:
            self.options_sampling.append((key,key))
        self.options_ori.clear()
        for key in self.orientations:
            self.options_ori.append((key,key))

    def _on_replot(self,b):
        # Plot experiment geometry 
        # (red arrow shows motion of sampling points in stationary sample)
        
        # read selected orientation
        self.on_change()
        orientation = self.orientations[self._widgets['ori'].value]
        sampling = self.sampling_list[self._widgets['sampling'].value]
        comm.set_sampling(sampling)
        comm.set_geometry(orientation)

       #  set_input(orientation, sampling)
        self._out.clear_output()
        with self._out:
            comm.plot_scene(self._widgets['nrec'].value, 
                            filename='', 
                            rang=2*[self._widgets['range'].value],
                            proj=self._widgets['proj'].value)    
        
    def update_widgets(self, data):
        self._update_options()
        self._widgets['sampling'].options = self.options_sampling
        self._widgets['ori'].options = self.options_ori
        for key in self._widgets:
            if key in data:
                self._widgets[key].value = data[key]

    def update_values(self):
        """Update input data from widgets.""" 
        for key in self._widgets:
            self._values[key] = self._widgets[key].value
        
    
            
    def show(self):
        lbl = create_header('Plot scene')
        box1 = ipy.HBox([self._widgets['sampling'], 
                         self._widgets['ori'], 
                         self._widgets['nrec']])
        btn_replot = ipy.Button(description='Replot',
                                layout=ipy.Layout(width='80px'))
        btn_replot.on_click(self._on_replot)
        box2 = ipy.HBox([self._widgets['proj'], 
                         self._widgets['range'], 
                         btn_replot])
        box = ipy.VBox([box1, box2])
        display(lbl)
        display(box)
        display(self._out)


class UI_workspace(UI_base):
    """
     Input for handling workspace directories.
     
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
                self._widgets[key] = PathInput(self._values, key, 
                                                  label=lbl, 
                                                  tooltip=ttp, 
                                                  output=self._out)
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
                wdg = self._widgets[key]
                wdg.set_value(data[key])  
        self.update_values()
        self.update_workspace()

    def update_values(self):
        """Update input data from widgets.
        
        NOTE: The widgets automatically update _values so this method is
        not in fact needed. It is provided only for consistency with the other 
        UI widtget collections.
        
        """ 
        for key in self._widgets.keys():
            self._values[key] = self._widgets[key].get_value()
    
    def show(self):
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
        box = ipy.VBox(lst)
        display(box)
        
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



 
 
#%% Top level UI class

class UI():
    def __init__(self):
        self.wk = dataio.workspace()
        self._out = ipy.Output(layout=ipy.Layout(width='100%'))
        # setup dialog
        self.setup = {}
        self.setup['workspace'] = UI_workspace()
        self.setup['shape'] = UI_shape(select=shapes.Tube)
        self.setup['geometry'] = UI_orientation()
        inpdir = self.setup['workspace'].get_values()['data']
        self.setup['sampling'] = UI_sampling()
        # initialize
        self.setup['geometry'].add_orientation('radial', update_table=True)
        sampl = {'path':inpdir, 'file':'events_S_1mm.dat', 'nev':3000}
        self.setup['sampling'].add_sampling('1x1x5mm',sampl, 
                                            update_table=True )
        # execution dialogs
        self.exec = {}
        self.exec['scene'] = UI_plot_scene(
            self.setup['sampling'].get_list(),
            self.setup['geometry'].get_list())
        self.exec['scene'].set_on_change(self._change)
    
    def display(self):
        self.setup['workspace'].show()
        self.setup['shape'].show()
        self.setup['geometry'].show()
        self.setup['sampling'].show()
        self.exec['scene'].show()
        display(self._out)
    
    def _change(self, obj, **kwargs):
        if isinstance(obj, UI_plot_scene):
            self._update_setup()
    
    def _update_setup(self):
        comp = self.setup['shape'].create_shape()
        comm.set_shape(comp)

    def get_values(self):
        out = {'setup':{}, 'exec':{}}
        for key in self.setup:
            out['setup'][key] = self.setup[key].get_values()
        for key in self.exec:
            out['exec'][key] = self.exec[key].get_values()
        return out

    def save(self, filename=''):
        out = {'ui': self.get_values} 
        txt = json.dumps(out,indent=4)
        if filename:
            wkp = self.wk.get_paths(keys=['work'])
            file = wkp['work'].joinpath(filename)
            f = open(file, 'w')
            f.write(txt)
            f.close()
        else:
            with self._out:
                print(txt)

            
        
