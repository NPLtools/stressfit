"""Widgets and classes built on ipywidgets.

Defines basic units for notebook interface.

Created on Tue Aug 15 13:44:06 2017
@author: Jan Saroun, saroun@ujf.cas.cz
"""
import ipywidgets as ipy
from traitlets import traitlets
from pathlib import Path as _Path
import numpy as np
import ipysheet
from IPython.display import display

def choose_path(initialdir=None):
    """Choose directory dialog using tkinter."""
    from tkinter import Tk, filedialog
    root = Tk()
    root.withdraw() # Hides small tkinter window.
    root.attributes('-topmost', True)
    open_file = filedialog.askdirectory(initialdir=initialdir)
    return open_file

def choose_file(initialdir=None, initialfile=None, filetypes=None, **kwargs):
    """Open file dialog using tkinter."""
    from tkinter import Tk, filedialog
    root = Tk()
    root.withdraw() # Hides small tkinter window.
    root.attributes('-topmost', True)
    ftall = ("All files","*.*")
    ftypes = []
    if filetypes is not None:
        ftypes.append(list(filetypes))
    ftypes.append(ftall)
    open_file = filedialog.askopenfilename(initialdir=initialdir, 
                                           initialfile=initialfile,
                                           filetypes=list(ftypes),
                                           **kwargs)
    return open_file


def choose_file_save(initialdir=None, initialfile=None, filetypes=None, **kwargs):
    """Open save_as file dialog using tkinter."""
    from tkinter import Tk, filedialog
    root = Tk()
    root.withdraw() # Hides small tkinter window.
    root.attributes('-topmost', True)
    ftall = ("All files","*.*")
    ftypes = []
    if filetypes is not None:
        ftypes.append(list(filetypes))
    ftypes.append(ftall)
    open_file = filedialog.asksaveasfilename(initialdir=initialdir, 
                                           initialfile=initialfile,
                                           filetypes=list(ftypes),
                                           **kwargs)
    return open_file


    
#%% Adapted widgets
# subclasses of ipywidgets: provide names to input widgets, and a button with value.
class SText(ipy.Text): 
    """Text with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SText, self).__init__(*args, **kwargs)
        self.name = name

class SBoundedFloatText(ipy.BoundedFloatText): 
    """BoundedFloatText with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SBoundedFloatText, self).__init__(*args, **kwargs)
        self.name = name

class SBoundedIntText(ipy.BoundedIntText): 
    """BoundedIntText with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SBoundedIntText, self).__init__(*args, **kwargs)
        self.name = name

class SDropdown(ipy.Dropdown): 
    """Dropdown with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SDropdown, self).__init__(*args, **kwargs)
        self.name = name

class SIntText(ipy.IntText): 
    """IntText with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SIntText, self).__init__(*args, **kwargs)
        self.name = name

class SFloatText(ipy.FloatText): 
    """FloatText with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SFloatText, self).__init__(*args, **kwargs)
        self.name = name


class SButton(ipy.Button):
    """Button with value trait."""
    
    def __init__(self, value='', *args, **kwargs):
        super(SButton, self).__init__(*args, **kwargs)
        self.add_traits(value=traitlets.Any(value))

class SCheckbox(ipy.Checkbox):
    """Check box with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SCheckbox, self).__init__(*args, **kwargs)
        self.name = name

class SRadioButtons(ipy.RadioButtons):
    """RadioButton box with name atribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SRadioButtons, self).__init__(*args, **kwargs)
        self.name = name


#%% Constructors of input elements
# Coponents for single value input with extended formatting features.

def create_input_float(name='', label='Float number',value=0.0, 
                       lmin=-1e30, lmax=1e30, step=0.1,
                       width_label=100, width_num=100, unit='px'):
    """Create named float input."""
    width = '{}{}'.format(width_label+width_num, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = SBoundedFloatText(name=name, description=label, value=value,
                                min=lmin, max=lmax, step=step, 
                                style=label_style,
                                layout=ipy.Layout(width=width))
    return inp

def create_input_int(name='', label='Integer',value=3000, 
                       lmin=None, lmax=None, step=1,
                       width_label=150, width_num=100, unit='px'):
    """Create named integer input."""
    width = '{}{}'.format(width_label+width_num, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = SBoundedIntText(name=name, description=label, value=value,
                                min=lmin, max=lmax, step=step, 
                                style=label_style,
                                layout=ipy.Layout(width=width))
    return inp


def create_header(text, color='#35446B', size='+1'):
    """Return header as HTML widget."""
    font_fmt = "<b><font color='{}' size='{}'>{}</font></b>"
    html = font_fmt.format(color,size,text)
    layout = ipy.Layout(padding='5px 20px 5px 0px')
    lbl = ipy.HTML(value=html, layout=layout)
    return lbl

   
def create_select(name='', options=[], label='', value=None, index=0, 
                  width_label=100, width_drop=60, unit='px'):
    """
    Create Dropdown widget for given options.

    Parameters
    ----------
    name: str
        Input variable name.
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
    sel = SDropdown(name=name, description=label,options=options,value=vini,
                           layout=ipy.Layout(max_width=max_width),
                           style=label_style)
    return sel

def create_checkbox(name='', label='Select',value=False, indent=False,
                    width_label=150, width_box=100, unit='px', tooltip=''):
    """Create named checkbox."""
    width = '{}{}'.format(width_label+width_box, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = SCheckbox(name=name, description=label, value=value, indent=indent,
                                style=label_style,
                                layout=ipy.Layout(width=width),
                                description_tooltip =tooltip)
    return inp

def create_text(name='', label='Text',value='', indent=False,
                width_label=150, width_box=100, unit='px', tooltip=''):
    """Create named text input."""
    width = '{}{}'.format(width_label+width_box, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = SText(name=name, description=label, value=value,
                continuous_update=False,
                style=label_style,
                layout=ipy.Layout(width=width),
                description_tooltip=tooltip)
    return inp

#%% Single-value input containers

class BasicInput():
    """Basic input class with name and callback on value change."""
    
    def __init__(self, name, value):
        self.name = name
        self._value = value
        self._observe = None
        self._call_enabled = True

    def _call_observe(self):
        """Invoke callback. To be used when value changes."""
        if self._observe and self._call_enabled:
            self._observe({'name':self.name, 'value': self.value})
    
    def observe(self, func):
        """Define callback function invoked when value changes."""
        self._observe = func
    
    @property
    def value(self):
        """Value getter."""
        return self._value
    
    @value.setter
    def value(self, value):
        """Value setter."""
        self._value = value


class DirInput(BasicInput):
    """Directory input as Text widget, open button and tooltip.
    
    Parameters
    ----------
    name: str
        Variable name
    path : str
        Path string
    label : str
        Label to appear before the text input.
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    
    Properties
    ----------
    value : :obj:`pathlib.Path`
        Directory path name.
    
    """
    
    def __init__(self, name='', path='', label='', tooltip='', 
                 width_label=150, width_button=50):
        super().__init__(name, path)
        self.label = label
        self.tooltip = tooltip
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = SText(value=self._value, name=name, 
                            layout=ipy.Layout(width='100%'),
                            continuous_update=False)
        self.btn = ipy.Button(description='', value=self.txt.value, 
                               layout=ipy.Layout(width=w_btn),
                               tooltip=self.tooltip, icon='folder-open')
        # add events
        self.btn.on_click(self._on_button)
        self.txt.observe(self._on_text,'value', type='change')


    def observe(self, func):
        """Define callback function invoked when value changes."""
        self._observe = func

    @property
    def value(self):
        """Value getter."""
        return self._value
    
    @value.setter
    def value(self, value):
        """Value setter."""
        # temporary switch off calling of observe.
        en = self._call_enabled
        self._call_enabled = False
        if isinstance(value, _Path):
            self.txt.value = value.as_posix()
        else:
            self.txt.value = str(value)
        self._call_enabled = en
            
    def _on_button(self,obj):
        """Open path dialog and update text input."""
        s = choose_path(initialdir=self.txt.value)
        if s:
            self.txt.value = s
    
    def _on_text(self, change):
        """Text change - update value."""
        if change['name']=='value':
            s = change['new']
            if s:
                self._value = s
                self._call_observe()

    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui
   

class FileInput(BasicInput):
    """File input as Text widget, open button and tooltip.
    
    Parameters
    ----------
    name: str
        Variable name        
    file : str
        Filename - full path or basename.
    path : str
        Optional directory name. 
    label : str
        Label to appear before the text input.
    fileonly : bool
        If true, value returns only filename as a string. 
        The path cannot be changed. Only files relative to the given path 
        are accepted. 
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    
    """
    
    def __init__(self,  name='', file='',path='', label='', fileonly=False,
                 tooltip='', width_label=150, width_button=50):
        
        super().__init__(name, '.')
        self._path = path
        self._fileonly = fileonly
        self._set_file(file)
        self.label = label
        self.tooltip = tooltip
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = ipy.Text(self._file,layout=ipy.Layout(width='100%'),
                            continuous_update=False)
        self.btn = ipy.Button(description='',  
                               layout=ipy.Layout(width=w_btn),
                               tooltip=self.tooltip, icon='folder-open')
        # add events
        self.btn.on_click(self._on_button)
        self.txt.observe(self._on_text,'value', type='change')


    @property
    def disabled(self):
        """Return disabled status."""
        return self.txt.disabled

    @disabled.setter
    def disabled(self, value):
        """Set disable status."""
        self.txt.disabled=value
        self.btn.disabled=value

    @property
    def value(self):
        """Return file name.
        
        If `self._fileonly`, then return just filename 
        as string. Otherwise return dict with `file` and `path` items. 
        
        """
        if self._fileonly:
            return self._value['file']
        else:
            return self._value
    
    @value.setter
    def value(self, value):
        """Set file and path names.
        
        Get file name from value and set new file and path names by
        callinf `self._set_file()`. If value contains `path`, update 
        `self._path` value first.
        
        Parameters
        ----------
        value : dict or str
            dict should contain `file` and `path` items
            str is iterpreted as file name. 
         
        """
        _path = self._path
        if isinstance(value, dict):
            if 'path' in value and value['path'] is not None:
                self._path = value['path']
            f = value['file']
        else:
            f = value
        try:    
            self._set_file(f)
        except Exception as e:
            self._path = _path
            raise e
        # temporary switch off calling of observe.
        en = self._call_enabled
        self._call_enabled = False
        self.txt.value = self._file
        self._call_enabled = en
    
    def _set_file(self, file):
        """Set new file name.
        
        1. Define full name. If file is relative, join it with 
           `self._path`.  
        2. Try to set `self._file` relative to `self._path`. If not possible, 
           then: if not `fileonly`,  define new `self._path`, else
           raise exception.
        """
        f = _Path(file)
        p = _Path(self._path)
        if f.is_absolute():
            fullname = f
        else:
            fullname = p.joinpath(file)
        try:
            _file = fullname.relative_to(p).as_posix()
        except:
            if not self._fileonly:
                p = f.parent
                _file = f.name
            else:
                msg = 'File name must be relative to {}'
                raise Exception(msg.format(p.as_posix()))
        self._fullname = fullname
        self._path = p.as_posix()
        self._file = _file
        self._value = {'path': self._path, 'file': self._file}                    
        
    def _on_button(self,ex):
        """Open path dialog and save result."""
        s = choose_file(initialdir=self._fullname.parent.as_posix(), 
                        initialfile=self._fullname.name)
        if s:
            try:
                self._set_file(s)
                self.txt.value = self._file
            except Exception as e:
                print(e)
             
    def _on_text(self, change):
        """Text change - update value."""
        if change['name']=='value':
            s = change['new']
            if s is not None:
                try:
                    self._set_file(s)
                    self._call_observe()
                except Exception as e:
                    print(e)
                                 
    def ui(self, width='100%', border='none'):
        """Return input widget as HBox with flex layout."""
        layout = ipy.Layout(display='flex', flex_flow='row', 
                                border=border, width=width)
        ui = ipy.HBox([self.lbl, self.txt, self.btn], layout=layout)
        return ui

    
class ArrayInput(BasicInput):
    """Array input set.
    
    Input set for a numerical array formatted as GridBox. Includes
    label, a row of numeric inputs and tooltip text.
    Scalar value is also accepted.
    
    Parameters
    ----------
    name: str
        Variable name  
    value: array_like
        Input array. The legth defines number of input widgets.
    label: str
        Label text on the left.
    hint: str
        Tooltip text on the right.
    step: float
        For float input,  defines step. 
    isInt: bool
        If true, the aray elements are treated as integers.
    width_label : int
        Width of label text in given units.
    width_num :int
        Width of numeric input widget in given units.
    unit : str
        Width units.
    """
    
    def __init__(self, name='', value=[0.0,0.0,0.0], label='array', 
                 hint='define an array', step=0.1, isInt=False,
                 width_label=150, width_num=80, unit='px',
                 lmin=None, lmax=None):
        super().__init__(name, value)
        if isinstance(value, (float,int)):
            self._value = [value]
            self._is_scalar = True
        else:
            self._value = value
            self._is_scalar = False
        self.label = label
        self.hint = hint
        self.isInt = isInt
        self.inputs = []
        if lmin is None:
            lmin=-1e30
        if lmax is None:
            lmax=1e30    
        x = [width_label,unit,len(self._value), width_num, unit]
        self._tpl = '{}{} repeat({},{}{}) auto'.format(*x)
        for j in range(len(self._value)):
            w = '{}{}'.format(width_num,unit)
            layout=ipy.Layout(width=w)
            if self.isInt:
                step = max(1,int(step))
                b = SBoundedIntText(name=j,description='',value=int(self._value[j]),
                                    step=step, layout=layout,
                                    min=lmin,max=lmax,
                                    continuous_update=False)
            else:
                b = SBoundedFloatText(name=j,description='',value=float(self._value[j]),
                                      step=step, layout=layout,
                                      min=lmin,max=lmax,
                                      continuous_update=False)
            self.inputs.append(b)
            b.observe(self._on_change,'value', type='change')

    @property
    def disabled(self):
        """Return disable status."""
        return self.inputs[0].disabled

    @disabled.setter
    def disabled(self, value):
        """Set disable status."""
        for b in self.inputs:
            b.disabled=value

    @property
    def value(self):
        """Get values from the widgets."""
        if self._is_scalar:
            arr = self._value[0]
        else:
            arr = self._value
        return arr 
    
    @value.setter
    def value(self, value):
        """Set values to the widgets."""
        if self._is_scalar:
            if not (isinstance(value, float) or isinstance(value, int)):
                raise Exception('Required scalar value, not {}'.format(value))
            self._value[0] = value
        else:
            ni = len(self.inputs)
            if len(value) != len(self.inputs):
                nd = len(value)
                raise Exception('Array lengths do not match: {},{}'.format(nd,ni))
            for i in range(ni):
                self._value[i] = value[i]
        en = self._call_enabled
        self._call_enabled = False
        self._update_widgets()
        self._call_enabled = en
                
    def _update_widgets(self):
        """Set values to the widgets."""
        ni = len(self.inputs)
        for i in range(ni):
            self.inputs[i].value = self._value[i]        
    
    def _on_change(self,change):
        """Value change - update value."""
        if change['name']=='value':
            s = change['new']
            if s is not None:
                idx = change['owner'].name 
                self._value[idx] = s
                self._call_observe()
    
    def ui(self):
        """Return the input container."""
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


class SelectInput(BasicInput):
    """Dropdown input set.
    
    Input set for an option list formatted as GridBox. Includes
    label, a Dropdown and tooltip text.
    
    Parameters
    ----------
    name: str
        Variable name 
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
    
    def __init__(self, name='',  options=[], label='', hint='', value=None, 
                 index=0, width_label=150, width_drop=60, unit='px'):
        super().__init__(name, value)
        self.options = options
        self.label = label
        self.hint = hint
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
                                      value=vini, layout=layout,
                                      continuous_update=False)
        self.input.observe(self._on_change,'value', type='change')
    
    @property
    def value(self):
        """Get selected value."""
        return self._value
    
    @value.setter
    def value(self, value):
        """Set value to the widget."""
        val = self._value
        self._value = value
        try:
            self._update_widgets()
        except  Exception as e:
            self._value = val
            raise e
        
    def set_index(self, index):
        """Set value by index to the widget."""
        if index >= len(self.options):
            msg = 'Index exceeds options length: {}>{}'
            raise Exception(msg.format(index,len(self.options)))
        try:
            self._call_enabled = False
            self.input.index = index 
        finally:
            self._call_enabled = True
        
    def _update_widgets(self):
        """Set values to the widgets."""
        en = self._call_enabled
        try:
            self._call_enabled = False
            self.input.value = self._value
            self._call_enabled = en
        except Exception as e:
            msg = 'Cannot set value to options: {}\n{}'
            self._call_enabled = en
            raise Exception(msg.format(self._value,e))

    def get_index(self):
        """Get selected index."""
        return self.input.index

    def _on_change(self,change):
        """Value change - update value."""
        if change['name']=='value':
           self._value = change['owner'].value
           self._call_observe()
                
    def ui(self):
        """Return the input container."""
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


class DistTable(BasicInput):
    """Table with xy nodes describing a free function.
    
    Beased in ipysheet, https://github.com/QuantStack/ipysheet/.
    
    Displays and handles a list of x,y coordinates and fix attributes. 
    It serves to describe a free distribution functions for fitting.
    
    To display the widget, call display(ui()), and then call redraw() 
    in order to display the table itself.
    
    Use the property 'value' to set or retrieve table data as dict. 
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        x, y, fix_x amd fix_y lists. 
        x,y are lists of float numbers. dix_x,y are boolean values.
    nx: int
        Number of nodes for default table.
    x_range: list(2)
        x-range for default table.
        
    """
    
    _headers = ['x','fix_x','y','fix_y']   
    def __init__(self, name='', value=None, nx=6, x_range=[0.0, 10.0],
                 border='1px solid', width='auto'):
        super().__init__(name, value)
        self.num_format = '0.0'
        self._cells = {}
        self._row_select = []      
        self._data = {}
        self._border = border
        self._width = width
        
        self._set_default_data(nx=nx, x_range=x_range)
        if isinstance(value, dict):
            self.import_data(value)
        self._data_orig = self._copy_data()
        layout = ipy.Layout(min_width=self._width, height='auto', 
                            border='none', 
                            margin='0px 0px 0px 0px')
        self._out = ipy.Output(layout=layout)
    
    @property
    def value(self):
        """Value getter."""
        return self.export_data()
    
    @value.setter
    def value(self, value):
        """Value setter."""
        self.import_data(value)
        self.redraw()
    
    def _set_default_data(self, nx=6, x_range=[0.0, 10.0]):
        [x1, x2] = x_range
        self._data['x'] = list(np.linspace(x1, x2, num=nx))
        self._data['fix_x'] = nx*[True]
        self._data['y'] = list(np.zeros(nx))
        self._data['fix_y'] = nx*[False]
        self._data_orig = self._copy_data()
    
    def _copy_data(self):
        """Return copy of table data."""
        data = {}
        for k in self._data:
            data[k] = self._data[k].copy()
        return data
    

    
    def _on_del_button(self, b):
        self._delete_rows()
        
    def _on_add_button(self, b):
        self._insert_rows()
    
    def _render_readonly(self, value):
        out = {}
        print(value)
        if value.read_only:
            out = {'backgroundColor' : "#EEEEEE"}
        return out
       
    def _create_sheet(self):
        """Create ipysheet from table data (self._data)."""
        nr = len(self._data['x'])  
        layout = ipy.Layout(min_width='auto', height='auto', 
                            border='none', 
                            margin='0px 0px 0px 0px')
        sheet = ipysheet.sheet(rows=nr, columns=len(DistTable._headers)+1,
                               row_headers = False,
                               column_headers = [' ']+DistTable._headers,
                               stretch_headers='none', layout=layout)
               
        stl = {'textAlign':'center'}
        gbcg = "#EEEEEE"
        
        # other data cells
        self._cells.clear()
        self._cells['x'] = ipysheet.column(1, self._data['x'][1:nr-1], 
                                                row_start=1, numeric_format=self.num_format)
        self._cells['fix_x'] = ipysheet.column(2, self._data['fix_x'][1:nr-1], 
                                                    row_start=1, type='checkbox', style=stl)
        self._cells['y'] = ipysheet.column(3, self._data['y'][1:nr-1], 
                                                row_start=1, numeric_format=self.num_format)
        self._cells['fix_y'] = ipysheet.column(4, self._data['fix_y'][1:nr-1], 
                                                    row_start=1, type='checkbox', style=stl)

        # x[0] and x[-1] must be always a fixed parameter        
        self._cells['00'] = ipysheet.cell(0,1, self._data['x'][0], type='numeric',
                                              numeric_format=self.num_format)
        ipysheet.cell(0,2, [None], read_only=True, background_color=gbcg);
        self._cells['02'] = ipysheet.cell(0,3, self._data['y'][0], type='numeric',
                                              numeric_format=self.num_format)
        self._cells['03'] = ipysheet.cell(0,4, self._data['fix_y'][0], type='checkbox', 
                                               style=stl, background_color="white")

        self._cells['10'] = ipysheet.cell(nr-1,1, self._data['x'][nr-1], type='numeric',
                                              numeric_format=self.num_format)
        ipysheet.cell(nr-1,2, [None], read_only=True, background_color=gbcg);
        self._cells['12'] = ipysheet.cell(nr-1,3, self._data['y'][nr-1], type='numeric',
                                              numeric_format=self.num_format)
        self._cells['13'] = ipysheet.cell(nr-1,4, self._data['fix_y'][nr-1], type='checkbox', 
                                               style=stl, background_color="white")
        
        # 1st column: check boxes for row selection
        ipysheet.cell(0, 0, [None], read_only=True, background_color=gbcg);
        self._row_select = ipysheet.column(0, (nr-2)*[False], row_start=1, 
                                           style=stl,
                                           background_color=gbcg)
        ipysheet.cell(nr-1, 0, [None], read_only=True, background_color=gbcg);
        self.sheet = sheet

    
    def _sheet_to_data(self):
        """Retrieve data from the sheet."""
        my_cells = ['00','02','03','10','12','13']
        if not all(k in self._cells.keys() for k in my_cells+DistTable._headers):
            raise Exception('DistTable: no sheet data.')
        data = {}
        data['x'] = [self._cells['00'].value] + self._cells['x'].value + [self._cells['10'].value]
        data['fix_x'] = [True] + self._cells['fix_x'].value + [True]
        data['y'] = [self._cells['02'].value] + self._cells['y'].value + [self._cells['12'].value]
        data['fix_y'] = [self._cells['03'].value] + self._cells['fix_y'].value + [self._cells['13'].value]
        self._data = data

    def _delete_rows(self):
        """Delete selected rows."""
        try:
            self._sheet_to_data()
            sel = self._row_select.value
            # collect data items to be deleted
            delc = {}
            for k in self._data:
                delc[k] = []
            for i in range(len(sel)):
                if sel[i]:
                    for k in self._data:
                        delc[k].append(self._data[k][i+1])
            # delete selected data items
            for k in delc:
                for d in delc[k]:
                    self._data[k].remove(d)
            self.redraw()
        except Exception as e:
            raise e

    def _insert_rows(self, before=False):
        """Insert rows after/before the selected ones."""
        try:
            self._sheet_to_data()
            sel = self._row_select.value
            d0 = self._copy_data()
            # collect selected indexes
            idx = []
            if before:
                di = 1
            else:
                di = 2
            for i in range(len(sel)):
                if sel[i]:
                    idx.append(i+di)            
            di = 0
            for i in idx:
                # insert values after i
                j = i + di
                x = 0.5*(d0['x'][i-1] + d0['x'][i])
                y = 0.5*(d0['y'][i-1] + d0['y'][i])
                self._data['x'].insert(j,x)
                self._data['y'].insert(j,y)
                self._data['fix_x'].insert(j,d0['fix_x'][i])
                self._data['fix_y'].insert(j,d0['fix_y'][i])
                di += 1
            self.redraw()
        except Exception as e:
            raise e

    def reset(self):
        """Reset to original value."""
        self.import_data(self._data_orig)
        self.redraw()
    
    def import_data(self, data, set_as_orig=False):
        """Import new table data."""
        if not all(k in data.keys() for k in DistTable._headers):
            raise Exception('DistTable import: missing keys.')
        nr = len(data['x'])
        for key in data:
            if nr != len(data[key]):
                raise Exception('DistTable import: unequal column lengths.')
            self._data[key] = data[key]
        if set_as_orig:
            self._data_orig = self._copy_data()
    
    def export_data(self):
        """Export sheet to data as dict."""
        self._sheet_to_data()
        return self._data
    
    def export_as_arrays(self):
        """Return data as 4 numpy arrays.
        
        Returns
        -------
        dict
            - x :  x-coordinates an array of float
            - y :  y-values an array of float
            - fix_x : fix flags for x values as an array of int
            - fix_y : fix flags for y values as an array of int
        """
        self._sheet_to_data()
        return {'x': np.array(self._data['x']), 
                'y': np.array(self._data['y']),
                'fix_x':np.array(np.array(self._data['fix_x']),dtype=int),
                'fix_y':np.array(np.array(self._data['fix_xy']),dtype=int)}
            
    def redraw(self):
        """Redraw the input table."""
        self._create_sheet()
        self._out.clear_output()
        with self._out:
            display(self.sheet)
    
    def ui(self):
        """Return the input container."""
        del_btn = ipy.Button(description='delete', 
                             layout=ipy.Layout(width='60px'), 
                             tooltip='delete row')
        del_btn.on_click(self._on_del_button)
        add_btn = ipy.Button(description='insert', 
                             layout=ipy.Layout(width='60px'), 
                             tooltip='insert row')
        add_btn.on_click(self._on_add_button)
        #lbl = ipy.Label('Rows: ')
        layout = ipy.Layout(width=self._width, border=self._border, 
                            margin='0px 0px 0px 0px')
        hb = ipy.HBox([del_btn, add_btn])
        self._create_sheet()
        out = ipy.VBox([hb, self._out],layout=layout)
        return out
      


class ScaleTable(BasicInput):
    """Table with distribution scale parameters.
    
    Beased in ipysheet, https://github.com/QuantStack/ipysheet/.
       
    Use the property 'value' to set or retrieve table data as dict. 
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        - keys : names of the values as list 
        - value : y-scale, y-offset and x-offset as list
        - fix : corresponding fix attributes as bool as list
    """
    
    _headers = ['keys','values','fix']
    
    def __init__(self, name='', value=None, fmt='0.00'):
        super().__init__(name, value)
        self.num_format = fmt
        self._cells = {}
        self._data = {}
        self._set_default_data()
        if isinstance(value, dict):
            self.import_data(value)
    
    @property
    def value(self):
        """Value getter."""
        return self.export_data()
    
    @value.setter
    def value(self, value):
        """Value setter."""
        self.import_data(value)
    
    def _set_default_data(self):
        self._data['keys'] = ['y-scale', 'y0', 'x0']
        self._data['values'] = [1.0, 0.0, 0.0]
        self._data['fix'] = [True, True, True]
        self._data_orig = self._copy_data()
    
        
    def _copy_data(self):
        """Return copy of table data."""
        data = {}
        for k in self._data:
            data[k] = self._data[k].copy()
        return data
    
                
    def _create_sheet(self):
        """Create ipysheet from table data (self._data)."""
        layout=ipy.Layout(width='auto', height='auto', border='none', 
                          margin='0px 0px 0px 0px')
        sheet = ipysheet.sheet(rows=3, columns=2,
                               row_headers = self._data['keys'],
                               column_headers = ['value', 'fixed'],
                               stretch_headers='none', layout=layout)
        stl = {'textAlign':'center'}
        vals = self._data['values']
        fixes = self._data['fix']
        self._cells.clear()
        self._cells['values'] = ipysheet.column(0, vals, row_start=0, 
                                               numeric_format=self.num_format)
        self._cells['fix'] = ipysheet.column(1, fixes, row_start=0, 
                                             type='checkbox', style=stl)
        self.sheet = sheet
    
    def _sheet_to_data(self):
        """Retrieve self._data from the sheet."""
        for key in self._cells:
            self._data[key] = self._cells[key].value
    
    def _data_to_sheet(self):
        """Set self._data to the sheet."""
        for key in self._cells:
            if key in self._data:
                self._cells[key].value = self._data[key]
        self.sheet.row_headers = self._data['keys']
        for c in self._cells:
            self._cells[c].send_state()

    def reset(self):
        """Reset to original value."""
        self.import_data(self._data_orig)
    
    def import_data(self,data, set_as_orig=False):
        """Import new data."""
        if not all(k in data.keys() for k in ScaleTable._headers):
            raise Exception('Invalid import data: missing keys.')
        nr = len(data['values'])
        for key in ScaleTable._headers:
            if nr != len(data[key]):
                raise Exception('Invalid import data: unequal array lengths.')
            self._data[key] = data[key]
        if set_as_orig:
            self._data_orig = self._copy_data()
        self._data_to_sheet()
    
    def export_data(self):
        """Export sheet to data as dict."""
        self._sheet_to_data()
        return self._data
    
    def export_as_arrays(self):
        """Return data as 2 numpy arrays.
        
        Returns
        -------
        dict
            - value : scale, y0 and x0 as an array of float
            - fix : fix flags as an array of int
        """
        self._sheet_to_data()
        vals = self._data['values']
        fixes = self._data['fix']
        return {'values': np.array(vals), 'fix':np.array(fixes,dtype=int)}
    
    def redraw(self):
        """Update table after manual change of parameters."""
        self._data_to_sheet()
        
    def ui(self):
        """Return the input widget."""
        self._create_sheet()
        return self.sheet
