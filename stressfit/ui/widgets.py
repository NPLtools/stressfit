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

# Information about working root directory is required to convert 
# relative paths to absolute and vice versa.


_workspace = {'root': _Path.cwd().as_posix()}

def set_workspace(workspace):
    """Set workspace information.

    Parameters
    ----------
    workspace : dict or str
        Information about user environment. 
        Currently only 'work' item with absolute path 
        to workspace root directory is accepted.
    """ 
    global _workspace
    if isinstance(workspace, dict) and 'work' in workspace:
        w = _Path(workspace['root'])
    elif isinstance(workspace, str):
        w = _Path(workspace)
    else:
        raise Exception('Only str and dict are accepted as workspace input')
    
    if w.is_absolute():
        _workspace['root'] = w
    else:
        raise Exception('Workspace directory must be an absolute path.')


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
    """Text with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SText, self).__init__(*args, **kwargs)
        self.name = name

class SBoundedFloatText(ipy.BoundedFloatText): 
    """BoundedFloatText with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SBoundedFloatText, self).__init__(*args, **kwargs)
        self.name = name

class SBoundedIntText(ipy.BoundedIntText): 
    """BoundedIntText with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SBoundedIntText, self).__init__(*args, **kwargs)
        self.name = name

class SDropdown(ipy.Dropdown): 
    """Dropdown with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SDropdown, self).__init__(*args, **kwargs)
        self.name = name

class SIntText(ipy.IntText): 
    """IntText with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SIntText, self).__init__(*args, **kwargs)
        self.name = name

class SFloatText(ipy.FloatText): 
    """FloatText with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SFloatText, self).__init__(*args, **kwargs)
        self.name = name


class SButton(ipy.Button):
    """Button with value trait."""
    
    def __init__(self, value='', *args, **kwargs):
        super(SButton, self).__init__(*args, **kwargs)
        self.add_traits(value=traitlets.Any(value))

class SCheckbox(ipy.Checkbox):
    """Check box with name attribute."""
    
    def __init__(self, name='', *args, **kwargs):
        super(SCheckbox, self).__init__(*args, **kwargs)
        self.name = name

class SRadioButtons(ipy.RadioButtons):
    """RadioButton box with name attribute."""
    
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
    """Basic input class with name and callback on value change.
    
    Do not use directly, implement own descendant classes.
    """
    
    class _WidgetLogger():
        def error(message):
            print(message)
        def exception(message):
            print(message)
        def info(message):
            print(message)
    
    def __init__(self, name, value, logger=None):
        self.name = name
        self._value = value
        self._observe = None
        self._call_enabled = True
        if logger==None:
            self._log = BasicInput._WidgetLogger()
        else:
            self._log = logger

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

    def ui(self):
        """Return the input widget container."""
        return None


class ComposedInput(BasicInput):
    """An input composed of multiple widgets for a set of values.
       
    Use the property 'value' to set or retrieve fit configuration data. 
    
    This is a base class which does nothing. Descendants should implement
    at least _create_values and _create_widgets methods.
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        Initial values.
    align : str
        vertical or horizontal alignment
    layout : obj
        Layout object to be passed to the flex container (either VBox or HBox)
    """
    
    def __init__(self, name='', value=None, align='vertical', 
                 layout = None, **kwargs):
        super().__init__(name, {}, **kwargs)
        self._widgets = {}
        self.align = align
        if layout:
            self.layout = layout
        else:
            self.layout = ipy.Layout()
        self._create_values()
        if isinstance(value, dict):
            self.value = value
    
    @property
    def value(self):
        """Value getter."""
        self._update_values()        
        return self._value
    
    @value.setter
    def value(self, value):
        """Value setter."""
        for key in value:
            if key in self._value:
                self._value[key] = value[key]
        self._update_widgets()
                
    def _create_values(self):
        """Define initial input data, one per widget.
        
        Called by constructor.
        """
        pass
    
    def _create_widgets(self):
        """Create widgets, one per value.
        
        Called by :meth:`ui` when creating the widgets. 
        """
        pass
    
    def _update_values(self):
        """Update values from widgets."""
        for key in self._widgets:
            if key in self._value:
                self._value[key] = self._widgets[key].value    
                
    
    def _update_widgets(self):
        """Update widgets from values."""
        for key in self._widgets:
            if key in self._value:
                self._widgets[key].value = self._value[key]

    def ui(self):
        """Return the input widgets."""
        self._create_widgets()
        wdg = []
        for key in self._widgets:
            if isinstance(self._widgets[key], BasicInput):
                wdg.append(self._widgets[key].ui())
            else:
                wdg.append(self._widgets[key])
        if self.align == 'vertical':
            out = ipy.VBox(wdg, layout=self.layout)       
        else:
            out = ipy.HBox(wdg, layout=self.layout)
        return out


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
                 width_label=150, width_button=50, **kwargs):
        super().__init__(name, path, **kwargs)
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

    
    @BasicInput.value.setter
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
     
        
    def fullname(self, path):
        """Return absolute path to the file."""
        p = _Path(path)
        if not p.is_absolute():
            p = _workspace['root'].joinpath(p)
        return p
        
    def _on_button(self,obj):
        """Open path dialog and update text input."""
        full = self.fullname(self.txt.value)
        s = choose_path(initialdir=full)
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
                 filetypes=None,
                 tooltip='', width='100%', width_label=150, width_button=50,
                 **kwargs):
        
        super().__init__(name, '.', **kwargs)
        self._path = path
        self._fileonly = fileonly
        self._filetypes = filetypes
        self._set_file(file, path=path)
        self.label = label
        self.tooltip = tooltip
        w_lbl = '{}px'.format(width_label)
        w_btn = '{}px'.format(width_button)
        self.lbl = ipy.Label(self.label, layout=ipy.Layout(min_width=w_lbl))
        self.txt = ipy.Text(self._file,layout=ipy.Layout(width=width),
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
            return self._file
        else:
            return self. _to_dict()
    
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
        try:
            if isinstance(value, dict):
                if 'path' in value and value['path'] is not None:
                    p = value['path']
                else:
                    p = self._path
                f = value['file']
                self._set_file(f, path=p)
            else:
                self._set_file(value)
        except Exception as e:
            self._log.exception(str(e))
        # temporary switch off calling of observe.
        en = self._call_enabled
        self._call_enabled = False
        self.txt.value = self._file
        self._call_enabled = en
    
    def fullname(self):
        """Return absolute path to the file."""
        p = _Path(self._path)
        if not p.is_absolute():
            p = _workspace['root'].joinpath(p)
        return p.joinpath(self._file)
    
    def _to_dict(self):
        """Return filename and path info as dict.
        
        Use current workspace info and try to derive relative paths to it.
        
        Set local _path and _file fields.
        
        Return
        ------
        dict
            path annd file values
        
        """       
        f = self.fullname()   
        try:
            rel = f.relative_to(_workspace['root'])
            path = rel.parent
            name = rel.name
        except:
            path = f.parent
            name = f.name
        # keep local variables updated, because workspace might have changed ... 
        self._path = path.as_posix()
        self._file = name
        return {'path': self._path, 'file': self._file}      
        
    def _set_file(self, file, path=None):
        """Set new file name.
        
        Derive also path using the path info and workspace root.

        """
        f = _Path(file)      
        self._file = f.name      
        if f.is_absolute():
            self._path = f.parent.as_posix()
        elif path:
            p = _Path(path).joinpath(f)
            self._path = p.parent.as_posix()        
        # _path is set as possibly absolute path
        # _to_dict() would tyr and replace it as relative to workspace root 
        out = self._to_dict()
        self._value = out
        
    def _on_button(self,ex):
        """Open path dialog and save result."""
        full = self.fullname()
        s = choose_file(initialdir=full.parent.as_posix(), 
                        initialfile=full.name,
                        filetypes=self._filetypes)
        if s:
            try:
                self._set_file(s)
                self.txt.value = self._file
            except Exception as e:
                self._log.exception(str(e))
             
    def _on_text(self, change):
        """Text change - update value."""
        if change['name']=='value':
            s = change['new']
            if s is not None and s != "":
                try:
                    self._set_file(s)
                    self._call_observe()
                except Exception as e:
                    self._log.exception(str(e))
                                 
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
    lmin : float
        Minimum value. None sets limit to 1e-30.
    lmax : float
        Maximum value. None sets limit to 1e30.
    """
    
    def __init__(self, name='', value=[0.0,0.0,0.0], label='array', 
                 hint='define an array', step=0.1, isInt=False,
                 width_label=150, width_num=80, unit='px',
                 lmin=None, lmax=None, **kwargs):
        super().__init__(name, value, **kwargs)
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
                self._log.error('Required scalar value, not {}'.format(value))
            self._value[0] = value
        else:
            ni = len(self.inputs)
            if len(value) != len(self.inputs):
                nd = len(value)
                self._log.error('Array lengths do not match: {},{}'.format(nd,ni))
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

class RangeInput(ArrayInput):
    """Set numeric range.
    
    Like ArratInput, but requires only 2-element value. 
    Provides range checking on values (value[0]<value[1]). 
    
    Parameters
    ----------
    name : str
        Name of this component. 
    value : list(2)
        Minimum and maximum values
    kwargs : dict
        Other named arguments passed to ArrayInput
    
    """
    
    def __init__(self, name='', value=[0.0, 1.0], **kwargs):         
        super().__init__(name=name, value=self._check_input(value), **kwargs)

    def _check_input(self, value):
        if not len(value)==2:
            raise Exception('RangeInput requires 2-element list as a value.')        
        if value[1]<value[0]:
            value[0], value[1] = value[1], value[0] 
        return value

    @property
    def value(self):
        """Get values from the widgets."""
        return self._value 
    
    @value.setter
    def value(self, value):
        """Set values to the widgets."""
        value = self._check_input(value)
        super(RangeInput, type(self)).value.fset(self, value)
        
    def _validate(self, change):
        iswp = [1, 0]
        sgn = [1, -1]
        i = change['owner'].name
        i2 = iswp[i]
        valid = sgn[i]*change['new'] < sgn[i]*self._value[i2]
        if not valid:
            change['owner'].value = change['old']
        return valid
    
    def _on_change(self, change):
        """Value change - update value."""
        if change['name']=='value':
            if self._validate(change):
                super()._on_change(change)

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
                 index=-1, width_label=150, width_drop=60, unit='px', **kwargs):
        super().__init__(name, value, **kwargs)
        self._options = options
        self.label = label
        self.hint = hint
        
        x = [width_label,unit, width_drop, unit]
        self._tpl = '{}{} {}{} auto'.format(*x)
        w = '{}{}'.format(width_drop,unit)
        layout = ipy.Layout(width=w)
        self.input = ipy.Dropdown(description='', options=self._options, 
                                      layout=layout,
                                      continuous_update=False)
        if value is not None:
            self.value = value
        else:
            self.index = index
        self.input.observe(self._on_change,'value', type='change')
        
    
    @property
    def options(self):
        """Get options."""
        return self._options
    
    @options.setter
    def options(self, value):
        self._options = value
        val = self._value # remember old value
        self.input.options = self._options
        # try to restore old value if possible
        if len(self._options)==0:
            return
        # get list of values from options
        if isinstance(self._options[0],str):
            values = list(self._options)
        else:
            values = list(dict(self._options).values())
        if val != self.input.value and val in values:
            self.input.value = val
        
        
    @BasicInput.value.setter
    def value(self, value):
        """Set value to the widget."""
        val = self._value
        self._value = value
        try:
            self._update_widgets()
        except  Exception as e:
            self._value = val
            self._log.exception(str(e))
       
    @property
    def index(self):
        """Return selected index."""
        return self.input.index
    
    @index.setter
    def index(self, index):
        """Set value by index to the widget."""
        if index>-1 and index<len(self._options):
            i = index  
        elif len(self._options)>0:
            i = 0
        else:
            i = -1
        if i>=0:
            try:
                self._call_enabled = False
                self.input.index = i 
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
            self._log.exception(msg.format(self._value,e))

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


class StatusBar(ComposedInput):
    """Display a horizontal bar with updatable values.
    
    Parameters
    ----------
    name : str
        Name of this component. 
    value : dict
        Any hash map with values. Must be provided at initiation. Later
        updates must include the same keys and value types.
    kwargs : dict
        Other named arguments passed to ComposedInput
    """
    
    def __init__(self, name='', value=None, **kwargs):
        if not isinstance(value, dict):
            raise Exception("StatusBar requires dict value on construction.")
        self._input = value    
        super().__init__(name, value, align='horizontal', **kwargs)
                
    def _create_values(self):
        self._value.clear()
        self._value.update(self._input)
 
    def _create_widgets(self):
        """Create value widgets."""
        n = len(self._value)
        keys = list(self._value.keys())
        wstr = '{:d}%'.format(int(80/n))
        for i in range(n):
            key = keys[i]
            wdg = ipy.Label(value='', layout=ipy.Layout(min_width=wstr)) 
         #                  layout=ipy.Layout(width=width))
            self._widgets[key] = wdg
        self._update_widgets()

    def _update_values(self):
        """Update values from widgets."""
        pass
    
    def _update_widgets(self):
        """Update widgets from values."""
        for key in self._widgets:
            if key in self._value:
                self._widgets[key].value = '{}: {:g}'.format(key, self._value[key])


class DistTable(BasicInput):
    """Table with xy nodes describing a free function.
    
    Based on ipysheet, https://github.com/QuantStack/ipysheet/.
    
    Displays and handles a list of x,y coordinates and fit attributes. 
    It serves to describe a free distribution functions for fitting.
    
    To display the widget, call display(ui()), and then call redraw() 
    in order to display the table itself.
    
    Use the property 'value' to set or retrieve table data as dict. 
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        x, y, fitx amd fity lists. 
        x,y are lists of float numbers. dix_x,y are boolean values.
    nx: int
        Number of nodes for default table.
    x_range: list(2)
        x-range for default table.
        
    """
    
    _headers = ['x','fitx','y','fity']   
    def __init__(self, name='', value=None, nx=6, x_range=[0.0, 10.0],
                 border='1px solid', label='distribution', **kwargs):
        super().__init__(name, value, **kwargs)
        self.num_format = '0.0'
        self._cells = {}
        self._row_select = []      
        self._data = {}
        self._border = border
        self._label = label
        
        self._set_default_data(nx=nx, x_range=x_range)
        if isinstance(value, dict):
            self.import_data(value)
        self._data_orig = self._copy_data()
        layout = ipy.Layout(width='100%', height='auto', border='none', 
                            padding='0px')
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
        self._data['fitx'] = nx*[1]
        self._data['y'] = list(np.zeros(nx))
        self._data['fity'] = nx*[0]
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
        if value.read_only:
            out = {'backgroundColor' : "#EEEEEE"}
        return out
       
    def _create_sheet(self):
        """Create ipysheet from table data (self._data)."""
        nr = len(self._data['x'])  
        layout = ipy.Layout(width='auto', height='auto', 
                            border='none', 
                            margin='0px', padding='0px')
        sheet = ipysheet.sheet(rows=nr, 
                               columns=len(DistTable._headers)+1,
                               row_headers = False,
                               column_headers = [' ']+DistTable._headers,
                               stretch_headers='none', layout=layout)
        sheet.column_width = [50] + len(DistTable._headers)*[50]
        stl = {'textAlign':'center'}
        gbcg = "#EEEEEE"
        
        # convert int to to bool
        fitx = list(map(bool,self._data['fitx']))
        fity = list(map(bool,self._data['fity']))
        
        # fill table cells, treat 1st and last columns separately
        self._cells.clear()
        
        self._cells['x'] = ipysheet.column(1, 
                                           self._data['x'][1:nr-1], 
                                           row_start=1, 
                                           numeric_format=self.num_format)
        self._cells['fitx'] = ipysheet.column(2, 
                                               fitx[1:nr-1], 
                                               row_start=1, 
                                               type='checkbox', 
                                               style=stl)
        self._cells['y'] = ipysheet.column(3, 
                                           self._data['y'][1:nr-1], 
                                           row_start=1, 
                                           numeric_format=self.num_format)
        self._cells['fity'] = ipysheet.column(4, 
                                               fity[1:nr-1], 
                                               row_start=1, 
                                               type='checkbox', 
                                               style=stl)

        # 1st row
        # 1st x must be always a fixed parameter        
        self._cells['00'] = ipysheet.cell(0, 1, 
                                          self._data['x'][0], 
                                          type='numeric',
                                          numeric_format=self.num_format)
        ipysheet.cell(0,2, [None], read_only=True, background_color=gbcg);        
        self._cells['02'] = ipysheet.cell(0, 3, 
                                          self._data['y'][0], 
                                          type='numeric',
                                          numeric_format=self.num_format)
        self._cells['03'] = ipysheet.cell(0, 4, 
                                          fity[0], 
                                          type='checkbox', 
                                          style=stl, 
                                          background_color="white")
        
        # last row
        # last x must be always a fixed parameter  
        self._cells['10'] = ipysheet.cell(nr-1, 1, 
                                          self._data['x'][nr-1], 
                                          type='numeric',
                                          numeric_format=self.num_format)
        ipysheet.cell(nr-1,2, [None], read_only=True, background_color=gbcg);
        self._cells['12'] = ipysheet.cell(nr-1, 3, 
                                          self._data['y'][nr-1], 
                                          type='numeric',
                                          numeric_format=self.num_format)
        self._cells['13'] = ipysheet.cell(nr-1, 4, 
                                          fity[nr-1], 
                                          type='checkbox', 
                                          style=stl, 
                                          background_color="white")
        
        # 1st column: check boxes for row selection
        # 1st and last row cannot be selected
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
            self._log.error('DistTable: no sheet data.')
        data = {}
        data['x'] = [self._cells['00'].value] + self._cells['x'].value + [self._cells['10'].value]
        data['y'] = [self._cells['02'].value] + self._cells['y'].value + [self._cells['12'].value]
        # convert bool to int
        fit = [False] + self._cells['fitx'].value + [False]   
        data['fitx'] = list(map(int,fit))
        fit = [self._cells['03'].value] + self._cells['fity'].value + [self._cells['13'].value]
        data['fity'] = list(map(int,fit))
        self._data = data

    def _delete_rows(self):
        """Delete selected rows."""
        try:
            self._sheet_to_data()
            sel = self._row_select.value
            # must have at least 3 nodes
            lensel = min(len(sel), len(self._data['x']) - 3)
            # collect data items to be deleted
            delc = {}
            for k in self._data:
                delc[k] = []
            for i in range(lensel):
                if sel[i]:
                    for k in self._data:
                        delc[k].append(self._data[k][i+1])
            # delete selected data items
            for k in delc:
                for d in delc[k]:
                    self._data[k].remove(d)
            self.redraw()
        except Exception as e:
            self._log.exception(str(e))

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
                self._data['fitx'].insert(j,d0['fitx'][i])
                self._data['fity'].insert(j,d0['fity'][i])
                di += 1
            self.redraw()
        except Exception as e:
            self._log.exception(str(e))

    def reset(self):
        """Reset to original value."""
        self.import_data(self._data_orig)
        self.redraw()
    
    def import_data(self, data, set_as_orig=False):
        """Import new table data."""
        if not all(k in data.keys() for k in DistTable._headers):
            self._log.error('DistTable import: missing keys.')
        nr = len(data['x'])
        for key in data:
            if nr != len(data[key]):
                self._log.error('DistTable import: unequal column lengths.')
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
            - fitx : fit flags for x values as an array of int
            - fity : fit flags for y values as an array of int
        """
        self._sheet_to_data()
        return {'x': np.array(self._data['x']), 
                'y': np.array(self._data['y']),
                'fitx':np.array(self._data['fitx'],dtype=int),
                'fity':np.array(self._data['fity'],dtype=int)}
            
    def redraw(self):
        """Redraw the input table."""
        with ipysheet.hold_cells():
            self._create_sheet()
        self._out.clear_output()
        with self._out:
            display(self.sheet)
    
    def ui(self):
        """Return the input container."""
        btn_padd = '0px, 5px, 0px, 5px'
        del_btn = ipy.Button(description='delete', 
                             layout=ipy.Layout(width='auto', padding=btn_padd), 
                             tooltip='delete row')
        del_btn.on_click(self._on_del_button)
        add_btn = ipy.Button(description='insert', 
                             layout=ipy.Layout(width='auto', padding=btn_padd), 
                             tooltip='insert row')
        add_btn.on_click(self._on_add_button)
        #lbl = ipy.Label('Rows: ')
        main_layout = ipy.Layout(display='flex',
                            width='28em', 
                            border=self._border, 
                            margin='0px')
        
        dist_layout = ipy.Layout(display='flex',
                            border='none', 
                            padding='0px')
        
        btn_layout = ipy.Layout(display='flex', 
                                padding='0.4em')

        
        self._create_sheet()
        
        btn_pos = 'left'
        if btn_pos =='left':
            btn_layout.width = '7em'
            hb = ipy.VBox([del_btn, add_btn], layout=btn_layout)
            dist = ipy.HBox([hb, self._out], layout=dist_layout)
        elif btn_pos =='right':
            btn_layout.width = '7em'
            hb = ipy.VBox([del_btn, add_btn], layout=btn_layout)
            dist = ipy.HBox([self._out, hb], layout=dist_layout)
        else:
            btn_layout.width = '15em'
            hb = ipy.HBox([del_btn, add_btn], layout=btn_layout)
            dist = ipy.VBox([hb,self._out], layout=dist_layout)            
        
        if self._label:
            lb = create_header(self._label, size='+0')
            out = ipy.VBox([lb, dist],layout=main_layout)
        else:
            out = ipy.VBox([dist], layout=main_layout)
        return out
      


class ScaleTable(BasicInput):
    """Table with distribution scale parameters.
    
    Based in ipysheet, https://github.com/QuantStack/ipysheet/.
       
    Use the property 'value' to set or retrieve table data as dict. 
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        - keys : names of the values as list 
        - value : y-scale, y-offset and x-offset as list
        - fit : corresponding fit attributes as bool as list
    """
    
    _headers = ['keys','values','fit']
    
    def __init__(self, name='', value=None, fmt='0.00', label='scale', 
                 **kwargs):
        super().__init__(name, value, **kwargs)
        self.num_format = fmt
        self._cells = {}
        self._data = {}
        self._set_default_data()
        self._label = label
        self.sheet = None
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
        self._data['fit'] = [1, 1, 1]
        self._data_orig = self._copy_data()
    
        
    def _copy_data(self):
        """Return copy of table data."""
        data = {}
        for k in self._data:
            data[k] = self._data[k].copy()
        return data
    
                
    def _create_sheet(self):
        """Create ipysheet from table data (self._data)."""
        layout=ipy.Layout(width='14em', height='15ex', 
                          border='none', 
                          margin='0px', padding='0px')
        sheet = ipysheet.sheet(key=self.name, rows=3, columns=2,
                               row_headers = self._data['keys'],
                               column_headers = ['value', 'fit'],
                               stretch_headers='all', layout=layout)
        #sheet.column_width = [70, 70, 70]
        stl = {'textAlign':'center'}
        vals = self._data['values']
        fixes = list(map(bool,self._data['fit']))
        
        self._cells.clear()
        self._cells['values'] = ipysheet.column(0, vals, row_start=0, 
                                               numeric_format=self.num_format)
        self._cells['fit'] = ipysheet.column(1, fixes, row_start=0, 
                                             type='checkbox', style=stl)
        self.sheet = sheet       
    
    def _sheet_to_data(self):
        """Retrieve self._data from the sheet."""
        for key in self._cells:
            if key=='fit':
                # convert bool to int ...
                self._data[key] = list(map(int,self._cells[key].value))
            else:
                self._data[key] = self._cells[key].value
                
    def _data_to_sheet(self):
        """Set self._data to the sheet."""
        if self.sheet:
            #with ipysheet.hold_cells():
            for key in self._cells:
                if key in self._data:
                    if key=='fit':
                        # convert int to bool ...
                        self._cells[key].value = list(map(bool,self._data[key]))
                    else:
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
            self._log.error('Invalid import data: missing keys.')
        nr = len(data['values'])
        for key in ScaleTable._headers:
            if nr != len(data[key]):
                self._log.error('Invalid import data: unequal array lengths.')
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
            - fit : fit flags as an array of int
        """
        self._sheet_to_data()
        vals = self._data['values']
        fixes = self._data['fit']
        return {'values': np.array(vals), 'fit':np.array(fixes,dtype=int)}
    
    def redraw(self):
        """Update table after manual change of parameters."""
        self._data_to_sheet()
        
    def ui(self):
        """Return the input widget."""
        self._create_sheet()
        dist_layout = ipy.Layout(display='flex',
                                 width='auto',
                                 border='none', 
                                 padding='0.4em',
                                 margin='0px 0px 0px 10px')
        if self._label:
            #lb = create_header(self._label, size='+0')
            lb = ipy.Label(value=self._label)
            out = ipy.VBox([lb, self.sheet], layout=dist_layout)
        else:
            out = ipy.VBox([self.sheet], layout=dist_layout)
        return out
        
        
        return self.sheet


class FitControl(ComposedInput):
    """Control fit settings.
       
    Use the property 'value' to set or retrieve fit configuration data. 
    
    Parameters
    ----------
    name: str
        Name of this component. 
    value : dict
        - maxiter : int, maximum iterations
        - guess : bool, run guest fit first 
        - loops :int,  maximum number of bootstrap loops
        - ar : float, regularization parameter    
    kwargs : dict
        Other named arguments passed to ComposedInput
    """
    
    def __init__(self, name='', value=None, **kwargs):
        super().__init__(name, value, align='vertical', **kwargs)
                
    def _create_values(self):
        self._value['maxiter'] = 200
        self._value['guess'] = True
        self._value['loops'] = 1
        self._value['ar'] = 3.0
        
    
    def _create_widgets(self):
        """Create value widgets."""
        maxiter = ArrayInput(name='maxiter', 
                             label='Iterations', 
                             hint='',
                             value=self._value['maxiter'],
                             isInt=True, step=10, lmin=0, lmax=1000,
                             width_label=80)
        loops = ArrayInput(name='loops', 
                           label='Loops', 
                           hint='Number of bootstrap loops',
                           value=self._value['loops'],
                           isInt=True, step=1, lmin=1, lmax=10,
                           width_label=80)
        ar = ArrayInput(name='ar', 
                        label='a_r', 
                        hint='Regularization parameter',
                        value=self._value['ar'],
                        width_label=80)     
        guess = create_checkbox(name='guess', 
                                label='Guess',
                                value=self._value['guess'], 
                                indent=False,
                                width_label=80, width_box=80,  
                                tooltip='Quick approximation before fit')
        self._widgets['maxiter'] = maxiter
        self._widgets['loops'] = loops
        self._widgets['ar'] = ar
        self._widgets['guess'] = guess


class StrainZeros(ComposedInput):
    """Control strain scaling.
    
    Controls setting of position encoder shift (z0) and zero strain (eps0).
       
    Use the property 'value' to set or retrieve fit configuration data. 
    
    Parameters
    ----------
    name : str
        Name of this component. 
    value : dict
        - z0 : float, encoder scale shift
        - eps0 : float,  zero strain in 1e-6
    kwargs : dict
        Other named arguments passed to ComposedInput
    """
    
    def __init__(self, name='', value=None, **kwargs):
        super().__init__(name, value, align='vertical', **kwargs)
                
    def _create_values(self):
        self._value['z0'] = 0.0
        self._value['eps0'] = 0.0
    
    def _create_widgets(self):
        """Create value widgets."""
        z0 = ArrayInput(name='z0', 
                             label=' ', 
                             hint='Encoder zero shift [mm]',
                             value=self._value['z0'],
                             width_label=80, logger=self._log)
        eps0 = ArrayInput(name='eps0', 
                           label=' ', 
                           hint='Zero strain [1e-6]',
                           value=self._value['eps0'],
                           width_label=80, logger=self._log)

        self._widgets['z0'] = z0
        self._widgets['eps0'] = eps0

class RegControl(ComposedInput):
    """Control of regularization loop.
    
    Set regularization parameters range and number of loops to scan. 
    
    Parameters
    ----------
    name : str
        Name of this component. 
    value : dict
        - range : float(2), range of regularization parameters
        - steps : int, number of steps
    kwargs : dict
        Other named arguments passed to ComposedInput
    """
    
    def __init__(self, name='', value=None, **kwargs):
        super().__init__(name, value, align='vertical', **kwargs)
                
    def _create_values(self):
        self._value['range'] = [0.0, 4.0]
        self._value['steps'] = 5
    
    def _create_widgets(self):
        """Create value widgets."""
        rang = RangeInput(name='range', 
                             label='Range', 
                             hint='Regularization range (log10 scale)',
                             value=self._value['range'],   
                             lmin=-10.0, lmax=20.0, step=0.5,
                             width_label=80, logger=self._log)
        steps = ArrayInput(name='steps', 
                           label='Steps', 
                           hint='Number of steps',
                           value=self._value['steps'],
                           isInt=True, step=1, lmin=2, lmax=64,
                           width_label=80, logger=self._log)

        self._widgets['range'] = rang
        self._widgets['steps'] = steps
                