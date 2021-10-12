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
    """Path or file input with text input, open button and tooltip.
    
    Saves result in the dict "result" under given key passed as arguments.
    
    For file input, set file=True. The `result` argument can be used to pass
    initial path and file strings under the keys `path` and `file`.
    
    Parameters
    ----------
    result: dict
        Dictionary where to store result.
    key: str
        Key where to store result.
    label: str
        Label to appear before the text input.
    tooltip: str
        A tooltip string (shows onthe button  mouse over event).
    file: bool
        If true, input should be a file name. 
    output: `ipywidgets.Output`
        Optional output widget is used to show error messages.
    
    """
    def __init__(self, result, key, label='', tooltip='', file=False, output=None):
        self.result = result
        self.key = key
        self.label = label
        self.tooltip = tooltip
        self.file = file
        self.output = output
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width='150px'))
        self.txt = ipy.Text(self.result[self.key], layout=ipy.Layout(width='100%'))
        self.btn = ValueButton(description='', value=self.txt.value, layout=ipy.Layout(width='10%'),
                             tooltip=self.tooltip, icon='folder-open')
        # link button and text
        self.dl = ipy.link((self.txt, 'value'), (self.btn, 'value'))
        # add events
        self.btn.on_click(self.on_button)
        self.txt.observe(self.on_text,'value', type='change')

    def on_button(self,ex):
        """Open file or path dialog and save result."""
        if self.file:
            f = None
            p = None
            if 'path' in self.result:
                p = self.result['path']
            if 'file' in self.result:
                f = self.result['file']
            s = choose_file(initialdir=p, initialfile=f)
        else:
            s = choose_path(initialdir=ex.value)
        if s:
            ex.value=s
            self.result[self.key] = s
            if self.output:
                self.output.clear_output()
    
    def on_text(self, change):
        """Text change - save result."""
        s = change['new']
        if s:
            self.result[self.key] = s
            if self.output:
                self.output.clear_output()
                        
    #def reset(self):
    #    self.txt.value = self.result[self.key]
        
    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui
    
    #def value(self):
    #    return self.b.value

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


class UI_shape(UI_base):
    """UI for sample shape definition.
    
    Parameters
    ----------
    select: str
        Selected shape, e.g. stressfit.shapes.Plate.
    
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
    """Input for a list of sample orientations.
    
    Each orientation is defined by four vectors:    
    
    angles : array_like(3)
        Euler YZY angles of sample orientation
    rotctr : array_like(3)
        Rotation centre in local sample coordinates [mm]
    scandir : array_like(3)
        Scan direction in local sample coordinates
    scanorig : array_like(3)
        Scan origin (position corresponding to zero scan position)
               
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
        if b.value in self._values:
            del self._values['list'][b.value]
            self.update_table()
    
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
        self._values['input'] = self._get_ori()

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
        

class UI_sampling(UI_base):
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
 
    
    def update_values(self):
        """Update input data from widgets.
        
        Does not change the list of defined sampling sets.

        """ 
        file = self._widgets['file'].get_value()
        self._values['input']['file'] = file['file']
        self._values['input']['path'] = file['path']
        self._values['input']['nev'] = self._widgets['nev'].value
       
    def update_widgets(self, data):
        """Update widgets from data.

        Parameters
        ----------
        data : dict
            Expected keys are [`file`,`path`,`nev`,`list`]
            
            `list` is a dict with meta-data [`file`,`path`,`nev`]
            for already loaded sampling sets. The method will try to
            reload these data and put them in self._sampling. 
            
            All keys are optional.

        """
        for key in self._keys:
            if key in data:
                self._values['input'][key] = data[key]
        if 'list' in data:
            lst = data['list']
            self._values['list'].clear()
            self._sampling.clear()
            for key in lst:
                self.add_sampling(key, lst[key])
            self._update_table()    
    
    def add_sampling(self, name, data, update_table=False):
        """Load sampling and add it to the list.
        
        Parameters
        ----------
        name : str
            Unique sampling name
        data : dict
            file, path and nev parameters
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
            self._widgets['file'].set_value(data)
            self._values['input'].update(data)
            self._values['list'][name]=data.copy()
            self._sampling[name]=sampling
            if update_table:
                self._update_table()         
    
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
        


class UI_plot_scene():
    def __init__(self, sampling_list, orientations, rang=16, proj=1, nev=3000):
        self.header = create_header('Plot scene')
        self.sampling_list = sampling_list
        self.orientations = orientations
        self.update_input()
        self.wdg_sampling = create_select(label='Select sampling', options=self.options_sampling, 
                                          width_label=100, width_drop=80)
        self.wdg_ori = create_select(label='Select orientation', options=self.options_ori, 
                                      width_label=130, width_drop=100)
        self.wdg_nrec = create_input_int(label='Events',
                               value=nev,
                               lmin=1000, lmax=100000, step=1000, 
                               width_label=150, width_num=100)
        #self.wdg_nrec = UI_sampling.events_input(width_label=130)
        self.wdg_projection = create_select(label='Projection: ', 
                                            options=[('z,y',0),('x,z',1),('x,y',2)],
                                            value = proj,
                                            width_label=100, width_drop=80)
        self.wdg_range = ipy.IntSlider(value=rang, min=3, max=100, step=1, description='Range: ')
        self.out = ipy.Output(layout=ipy.Layout(width='100%', border='none')) 
        
    def update_input(self):
        self.options_sampling = []
        #for i in range(len(self.sampling_list)):
        for key in self.sampling_list:
            self.options_sampling.append((key,key))
        self.options_ori = []
        for key in self.orientations:
            self.options_ori.append((key,key))

    def update_widgets(self):
        self.update_input()
        if self.wdg_sampling:
            self.wdg_sampling.options = self.options_sampling
        if self.sel_ori.options:
            self.wdg_ori.options = self.options_ori

    def on_replot(self,b):
        # Plot experiment geometry 
        # (red arrow shows motion of sampling points in stationary sample)
        
        # read selected orientation
        orientation = self.orientations[self.wdg_ori.value]
        sampling = self.sampling_list[self.wdg_sampling.value]
        #TODO

       #  set_input(orientation, sampling)
        self.out.clear_output()
        with self.out:
            comm.plot_scene(self.wdg_nrec.value, 
                            filename='', 
                            rang=[self.wdg_range.value,self.wdg_range.value],
                            proj=self.wdg_projection.value)        
            
    def show(self):
        lbl = create_header('Plot scene')
        box1 = ipy.HBox([self.wdg_sampling, self.wdg_ori, self.wdg_nrec])
        btn_replot = ipy.Button(description='Replot',layout=ipy.Layout(width='80px'))
        btn_replot.on_click(self.on_replot)
        box2 = ipy.HBox([self.wdg_projection, self.wdg_range, btn_replot])
        box = ipy.VBox([box1, box2])
        display(lbl)
        display(box)
        display(self.out)


class UI_workspace():
    
    def __init__(self, env):
        self.env = env
        self.out = ipy.Output(layout=ipy.Layout(width='100%'))
        self.path_input = PathInput(env, 'data', label='Input data', 
                                    tooltip='Path to input data.', 
                                    output=self.out)
        self.path_tables = PathInput(env, 'tables', label='Tables', 
                                     tooltip='Path to material tables etc.', 
                                     output=self.out)
        self.path_output = PathInput(env, 'output', label='Output', 
                                     tooltip='Path for data output.', 
                                     output=self.out)

       
    def validate_workspace(self,b):
        self.out.clear_output()
        with self.out:
            try:
                dataio.set_path(**self.env)
                dataio.validate_paths()
                print('OK')
            except Exception as e:
                print(e)
                
    def reset_workspace(self,b):
        self.out.clear_output()
        with self.out:
            dataio.reset_paths()
            self.env.update(dataio.get_paths())
            self.path_input.reset()
            self.path_tables.reset()
            self.path_output.reset()
            
    def show(self):
        appl = ipy.Button(description='Apply',layout=ipy.Layout(width='10%'))
        appl.on_click(self.validate_workspace)
        rst = ipy.Button(description='Reset',layout=ipy.Layout(width='10%'))
        rst.on_click(self.reset_workspace)
        hdr = create_header('Set workspace')
        
        layout = ipy.Layout(display='flex',flex_flow='row',
                                border='none',width='100%')
        appl_box = ipy.HBox([appl, self.out, rst], layout=layout)
        
        
        box = ipy.VBox([hdr, self.path_input.ui(), self.path_tables.ui(),
                            self.path_output.ui(), appl_box])
        display(box)
 
 
#%% Top level UI class

class UI():
    def __init__(self):
        self.input = {}
        self.input['env'] = dataio.get_paths()
        self.input['shapes'] = dataio.load_config('shapes')

    def set_input(orientation, sampling):
        """Apply user input to stressfit."""
       # TODO
        # comp = ui_shape.create_shape()
       # comm.set_shape(comp)
        keys =  ['scandir','scanorig','angles','rotctr']
        scan = {key: orientation[key] for key in keys}
        scan['sampling'] = sampling
        comm.set_scan(scan)

