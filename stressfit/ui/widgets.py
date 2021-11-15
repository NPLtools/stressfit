"""Widgets and classes built on ipywidgets.

Defines basic units for notebook interface.

Created on Tue Aug 15 13:44:06 2017
@author: Jan Saroun, saroun@ujf.cas.cz
"""
import ipywidgets as ipy
from traitlets import traitlets
from pathlib import Path as _Path

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
    lbl = ipy.HTML(value=html)
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
                                tooltip=tooltip)
    return inp

def create_text(name='', label='Text',value='', indent=False,
                width_label=150, width_box=100, unit='px', tooltip=''):
    """Create named text input."""
    width = '{}{}'.format(width_label+width_box, unit)
    label_style = {'description_width': '{}{}'.format(width_label,unit)} 
    inp = SText(name=name, description=label, value=value,
                                style=label_style,
                                layout=ipy.Layout(width=width),
                                tooltip=tooltip)
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
        self._call_enabled = False
        if isinstance(value, _Path):
            self.txt.value = value.as_posix()
        else:
            self.txt.value = str(value)
        self._call_enabled = True
            
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
    tooltip : str
        A tooltip string (shows onthe button  mouse over event).
    output : `ipywidgets.Output`
        Optional output widget is used to show error messages.
    width_label : int
        With of the label in px
    width_button : int
        Width of the load button in px.
    
    """
    
    def __init__(self,  name='', file='',path='', label='', tooltip='', 
                 output=None,
                 width_label=150, width_button=50):
        
        super().__init__(name, '.')
        self._path = path
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
        """Return file name as dict with file and path items."""
        return self._value
    
    @value.setter
    def value(self, value):
        """Set file name from dict with file and path items.
        
        The path item is optional. 
        """
        if 'path' in value and value['path'] is not None:
            self._path = value['path']
        self._set_file(value['file'])
        self._call_enabled = False
        self.txt.value = self._file
        self._call_enabled = True
    
    def _set_file(self, file):
        """Set new file name.
        
        If file is relative, join it with path info. 
        
        If file is absolute, update the path info.                
        """
        f = _Path(file)
        if f.is_absolute():
            p = f.parent
            self._fullname = f
            self._file = f.name
            self._path = p.as_posix()
        else:
            p = _Path(self._path)
            self._fullname = p.joinpath(file)
            try:
                self._file = self.fullname.relative_to(p).as_posix()
            except:
                self._file = file # this should not happen
        self._value = {'path': self._path, 'file': self._file}                    
        
    def _on_button(self,ex):
        """Open path dialog and save result."""
        s = choose_file(initialdir=self._fullname.parent.as_posix(), 
                        initialfile=self._fullname.name)
        if s:
            self._set_file(s)
            self.txt.value = self._file
            
    def _on_text(self, change):
        """Text change - update value."""
        if change['name']=='value':
            s = change['new']
            if s:
                self._set_file(s)
                self._call_observe()
                 
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
        self._call_enabled = False
        self._update_widgets()
        self._call_enabled = True
                
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
        try:
            self._call_enabled = False
            self.input.value = self._value
        except Exception as e:
            msg = 'Cannot set value to options: {}\n{}'
            self._call_enabled = True
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
