# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""
File input/output functions for data used by StressFit.

File input functions
--------------------
Use the `load_...` functions to load data in various formats.



By default, input data are searched for in package resources. There are three
types of such resource paths:
    
    `data`
        for input data (e.g. experimental data)
    `tables`
        for tables supplied with the package, such as lookup tables
        for material properties etc.
    `instruments`
        for instrument configurations.

Hints:
    - To define which of the resource paths to search in, use the ``kind`` parameter
      of the `load` functions. 
    - To define different default search paths in user space, use :func:`~.set_path`.
    - To override the default search paths, either provide file names as absolute 
      paths, or use the ``path`` parameter of the `load` functions.


File output functions
---------------------
Use the `save_...` functions to save data.

Hints:
    - By default, output data are saved in the current directory. 
    - To define another default output path, use :func:`~.set_path` with 
      the `output` parameter.
    - To override the default output paths, either provide file name as an absolute 
      path, or use the ``path`` parameter of the `save` functions.

Data formats
------------

This module saves/loads data as text files. Internally, the data can 
be represented by:

    :obj:`numpy.ndarry`
        A 2D array of float data. Use :func:`~.save_data` and 
        :func:`~.load_data` for file I/O.
    :class:`~.Table`
        A class defined by this module, which represents 2D tables as sets
        of rows. Each row is a list of objects. The class also provides 
        access to columns and meta data such as column/row headers and comment
        lines. 
        Use :func:`~.save_table` and :func:`~.load_table` for Table I/O.

    :obj:`dict` of `~.Param`
        A dictionary of `~.Param` objects as values. 
        Use :func:`~.save_params` and :func:`~.load_params` for Table I/O.
 
"""
import numpy as np
import datetime
import json
import traceback
try:
    from pygments import formatters, highlight, lexers
    _has_pygments = True
except ImportError:
    _has_pygments = False
try:
    import ipywidgets as ipy
    _has_ipy = True
except ImportError:
    _has_ipy = False
try:
    from IPython.display import display
    _has_display = True
except:
    _has_display = False
from pathlib import Path as _Path
import os
from functools import wraps
import copy
import logging
from colorama import Fore

__work = None
__setup = None
__log = None

# allowed workspace directories
_wks_keys = ['work', 'data', 'tables', 'instruments', 'output']
# workspace directories without root
_path_keys = _wks_keys.copy()
_path_keys.remove('work')


def logger():
    """Access Stressfit loggers."""
    global __log
    if __log is  None:
        __log = StressfitLogger()
    return __log


def _print_loggers(full=False):
    for k,v in  logging.Logger.manager.loggerDict.items()  :
        if full or k.strip().startswith('stressfit'):
            print('+ [%s] {%s} ' % (str.ljust( k, 20)  , str(v.__class__)[8:-2]) ) 
            if not isinstance(v, logging.PlaceHolder):
                for h in v.handlers:
                    print('     +++',str(h.__class__)[8:-2] )


class StressfitLogger():
    """Encapsulate python logging logger for Stressfit."""
    
# TODO progress logger
    def __init__(self):
                
        # define optional widget outputs for loggers
        self._msg = None # basic info, warning and error 
        self._exc = None # exceptions
        self._prog = None # progress
        self._short = None # short messages, 2-line output area 
        
        # message logger: info, warning, error
        self._lm = logging.getLogger('stressfit')
        self._lm.setLevel(logging.INFO)
        self._hm = _log_handler()
        self._hm.setFormatter(_log_formatter())
        self._lm.addHandler(self._hm)
        # optional file handler
        self._hm_file = None
        
        # exception logger
        self._lex = logging.getLogger('stressfit.exception')
        self._lex.setLevel(logging.ERROR)
        self._lex.propagate = False
        self._hex = _log_handler_exc()
        self._lex.addHandler(self._hex)
        # optional file handler
        self._hex_file = None
        
        # progress logger
        self._lprog = logging.getLogger('stressfit.progress')
        self._lprog.setLevel(logging.INFO)
        self._lprog.propagate = False
        self._hprog = _log_handler_prog()
        self._lprog.addHandler(self._hprog)
        # optional file handler
        self._hprog_file = None
        

    @property
    def output_msg(self): 
        """Output ipywidget, messages"""
        return self._msg
    
    @output_msg.setter
    def output_msg(self, value):
        self._msg = value
        if self._hm:
            self._hm.out = value

    @property
    def output_exc(self): 
        """Output ipywidget, exceptions"""
        return self._exc
    
    @output_exc.setter
    def output_exc(self, value):
        self._exc = value
        if self._hex:
            self._hex.output = value
        
    @property
    def output_prog(self): 
        """Output ipywidget, progress"""
        return self._prog
    
    @output_prog.setter
    def output_prog(self, value):
        self._prog = value
        if self._hprog:
            self._hprog.output = value

    @property
    def output_short(self): 
        """Output ipywidget, short messages"""
        return self._short
    
    @output_short.setter
    def output_short(self, value):
        self._short = value              

    @property
    def file_log(self):
        return self._hm_file is not None
    
    @file_log.setter
    def file_log(self, value:bool):
        if value and self._hm_file is None:
            self._hm_file = logging.FileHandler('messages.txt', mode='w')
            self._lm.addHandler(self._hm_file)
        elif not value:
            try:
                if self._hm_file:
                    self._hm_file.close()
                    if self._hm_file in self._lm.handlers:
                        self._lm.removeHandler(self._hm_file)
            finally:
                self._hm_file = None                
        return self._hm_file is not None
    
    
    def _print_short(self, message, level):
        """Print 1st line on the short output if defined."""
        if self._short and isinstance(message, str):
            # always clear it
            self._short.clear_output()
            lst = message.strip().split('\n')
            msg = lst[0]
            if level>=logging.ERROR:
                fore = Fore.RED
            elif level>=logging.WARNING:
                fore = Fore.MAGENTA
            else:
                fore = Fore.RESET
            with self._short:
                print(fore + msg + Fore.RESET)
    
    def info(self, message):
        """Print info message."""        
        if isinstance(message, list):
            message='\n'.join(message)
        self. _print_short(message, logging.INFO)
        self._lm.info(message)
        
    def warning(self, message):
        """Print warning."""   
        if isinstance(message, list):
            message='\n'.join(message)
        self. _print_short(message, logging.WARNING)
        self._lm.warning(message)
        
    def error(self, message):
        """Print error message."""
        if isinstance(message, list):
            message='\n'.join(message)
        self. _print_short(message, logging.ERROR)
        self._lm.error(message)
        
    def exception(self, message):
        """Print error message."""
        if isinstance(message, list):
            message='\n'.join(message)
        self. _print_short(message, logging.ERROR)
        self._lex.error(message)

    def progress(self, message):
        """Print error message."""
        if isinstance(message, list):
            message='\n'.join(message)
        self._lprog.info(message)
        
    def setLevel(self, level):
        """Set main logger level."""
        self._lm.setLevel(level)

    def add_handler(self, hnd):
        """Add another handler to the message logger."""
        self._lm.addHandler(hnd)  
    
    def remove_handler(self, hnd=None):
        """Remove given  handler from the message logger."""
        if hnd == None:
            hnd = self._hm
        self._lm.removeHandler(hnd)
    
    def remove_all_handlers(self):
        """Remove all handlers from gthe message logger."""
        while self._lm.hasHandlers(): 
            self._lm.removeHandler(self._lm.handlers[0])
    
    def reset_handlers(self):
        """Remove all handlers from gthe message logger except the default one."""
        self.remove_all_handlers()
        self._lm.addHandler(self._hm)         
        
    def clear(self, what='all'):
        """Clear output if defined.
        
        Parameters
        ----------
        what : str
            Output to be cleared:
                
            error 
                error output
            info 
                short info messages                
            msg
                list of messages
            prog
                progress info
            all
                all outputs
        """
        if what=='error' and self.output_exc:
            self.output_exc.clear_output()
        elif what=='prog' and self.output_prog:
            self.output_prog.clear_output()
        elif what=='info' and self.output_short:
            self.output_short.clear_output()
        elif what=='msg' and self.output_msg:
            self.output_msg.outputs = ()
            self.output_msg.clear_output()
        elif what=='all':
            if self.output_exc:
                self.output_exc.clear_output()
            if self.output_prog:
                self.output_prog.clear_output()
            if self.output_short:
                self.output_short.clear_output()
            if self.output_msg:
                self.output_msg.outputs = ()
                self.output_msg.clear_output()
            
                   
class _log_formatter(logging.Formatter):    
    FORMATS = {logging.WARNING: '%(levelname)s: %(message)s',
               logging.ERROR: '%(levelname)s: %(message)s',
               logging.CRITICAL: '%(levelname)s: %(message)s'}
    def format(self, record):
        if record.levelno in _log_formatter.FORMATS:
            log_fmt = _log_formatter.FORMATS.get(record.levelno)
        else:
            log_fmt = '%(message)s'
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
                    
class _log_handler(logging.Handler):
    """Class to redistribute python logging data."""
    
    def __init__(self, *args, **kwargs):
         # Initialize the Handler
         logging.Handler.__init__(self, *args)
         self.out = None
         # optional take format
         for key, value in kwargs.items():
             if "{}".format(key) == "output":
                 self.out = value

    def emit(self, record):
        """Overload of logging.Handler method."""
        # set collor for the level
        if record.levelno>=logging.ERROR:
            fore = Fore.RED
        elif record.levelno>=logging.WARNING:
            fore = Fore.MAGENTA
        else:
            fore = Fore.RESET

        # format record
        msg = self.format(record)
        if self.out:
            new_output = {
                'name': 'stdout',
                'output_type': 'stream',
                'text': fore+msg+Fore.RESET+'\n'
            }
            # add on top of the outputs
            self.out.outputs = (new_output, ) + self.out.outputs
        else:
            print(fore + msg + Fore.RESET)

class _log_handler_exc(logging.StreamHandler):
    """Exception logging handler."""
    
    def __init__(self, *args, **kwargs):
         # Initialize the Handler
         super().__init__(self, *args)
         self.output = None
         fmt = logging.Formatter('EXCEPTION: %(message)s')
         self.setFormatter(fmt)        

    def error_trace(self):
        """Format and print traceback from the last exception.""" 
        tb_text = "".join(traceback.format_exc()) 
        if _has_pygments:
            if _has_ipy and _has_display and self.output:
                lexer = lexers.get_lexer_by_name("pytb", stripall=True)
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
                with self.output:
                    display(t)
            else:
                # self.message('terminal16m')
                formatter = formatters.get_formatter_by_name('terminal16m')
                tb_colored = highlight(tb_text, lexer, formatter)
                print(tb_colored)
        elif self.output:

            with self.output:
                print(tb_text)
        else:
            print(tb_text)

    def emit(self, record):
        """Print to output if available."""
        msg = self.format(record)
        if self.output:
            self.output.clear_output()
            with self.output:
                print(Fore.RED + msg + Fore.RESET)
        else:
            print(Fore.RED + msg + Fore.RESET)
        self.error_trace()

class _log_handler_prog(logging.StreamHandler):
    """Progress logging handler."""
    
    def __init__(self, *args, **kwargs):
         # Initialize the Handler
         super().__init__(self, *args)
         self.output = None
         fmt = logging.Formatter('%(message)s')
         self.setFormatter(fmt) 

    def emit(self, record):
        """Print to output if available."""
        msg = self.format(record)
        if self.output:
            with self.output:
                print(msg )
        else:
            print(msg)

def workspace():
    """Access Workspace object.
    
    Return
    ------
    instance of :class:`dataio.Workspace`
        Workspace object which defines and maintains user 
        workspace directories.
    """
    global __work
    if __work is None:
        __work = Workspace()
    return __work


def _setup():
    """Access the instance of _Setup.
    
    _Setup handles application global settings and default values. It should
    remain "private" as there is nothing to be done on it by users.
    
    Return
    ------
    Instance of :class:`dataio._Setup`
    """
    global __setup
    if __setup is None:
        __setup = _Setup()
    return __setup




def _mk_path_relative(root, path):
    """Try to make given path relative to the root path.
    
    If already relative or root is not its parent, return the path unchanged.
    
    Parameters
    ----------
    root : pathlib.Path
        Parent path
    path : pathlib.Path
        Path to be converted.
        
    Return
    ------
    pathlib.Path:
        If possible, path relative to root, or unchanged path. 
    """
    out = path
    if path.is_absolute():
        try:
            out = path.relative_to(root)
        except:
            # print('cannot set relative:\n{}\n{}'.format(root,path))
            out = path
    return out


def _workspace_to_str(wks, absolute=False, keys=None):
    """Return workspace paths as dict of strings.
    
    Parameters
    ----------
    wks : dict
        Workspace paths.
    absolute : bool
        Return as absolute paths.
    keys : list
        If defined, return only listed paths. Otherwise
        return all workspace paths.
    
    Return
    ------
    dict
        Keys are the path types defined for the workspace.
        The values are path names returned by :func:`pahlib.Path.as_posix`.
    """
    env = {}
    if keys is None:
        keys = _wks_keys
    for key in keys:
        try:
            if absolute:
                if wks[key].is_absolute():
                    p = wks[key]
                else:
                    p = wks['work'].joinpath(wks[key])
            else:
                p = wks[key]
            env[key] = p.as_posix()
        except KeyError:
            print('Missing key in workspace data: {}'.format(key))
    return env

def _str_to_workspace(source):
    """Convert dict of strings to a dict of pathlib.Path.
    
    Makes paths relative to workspace root if possible.
    """
    out = {}
    for key in _wks_keys:
        try:
            out[key] = _Path(source[key])
        except KeyError:
            print('Missing required path in workspace: {}'.format(key))
    for key in _path_keys:
        out[key] = _mk_path_relative(out['work'], out[key])
    return out
            
def _create_workspace_tree(wks):
    """Create workspace directory tree.
    
    Creates directories only if they are childs of the workspace path.
    
    Parameters
    ----------
    wks : dict
        Workspace paths info.
    """
    for key in _path_keys:
        f = wks[key]
        p = _mk_path_relative(wks['work'], wks[key])
        # create only subdirectories
        if not p.is_absolute():
            f = wks['work'].joinpath(p)
            if not f.exists():
                f.mkdir(parents=True)
            elif f.is_file():
                msg = "File {} already exists.\nCannot create directory."
                raise Exception(msg.format(f.as_posix()))



class _Setup():
    """
    Class maintaining application settings and default values.
    
    Maintains a dict with version info, default paths, list of workspaces, etc.
    Provides load/save functions.
    
    A singleton instance of _Setup is cretaed internally when the function 
    _setup() is called for the first time.
    
    """
    
    _setup_dir = '.stressfit' # directory with workspace info
    _setup_file = 'setup.json' # file with workspace info
    
    def __init__(self):
        self._data = {}
        self._create()
    
    @property
    def content(self):
        return self._data
    
    @property
    def workspace(self):
        return self._data['workspace']
    
    def _create(self):
        """Create global configuration directories and files if needed."""
        # find application data folder
        ep = os.getenv('APPDATA')
        if ep:
            epath = _Path(ep).joinpath(_Setup._setup_dir)
        else:
            epath = _Path.home().joinpath(_Setup._setup_dir)
        self._appdir = epath
        # file with application setup info
        self._appsetup = epath.joinpath(_Setup._setup_file)
        # make sure the setup path exist
        if not self._appdir.exists():
            self._appdir.mkdir()
        elif self._appdir.is_file():
            self._appdir.unlink()
            self._appdir.mkdir()
        
        # is there already the setup file?
        if self._appsetup.exists():
            # load setup file
            with open(self._appsetup, 'r') as f:
                lines = f.readlines() 
            content = json.loads('\n'.join(lines))            
            self.verify(content)
            wks = _str_to_workspace(content['workspace'])
            content.update({'workspace':wks})
            self._data.update(content)
        else:
            # create content and save it
            cont = load_config('setup.json') # this is a template in resources
            self.verify(cont)
            self._data.update(cont)
            self.reset_paths()  # set paths to default (using resources)
            self.save()
    
    def verify(self, cont):
        """Verify setup file content."""
        out = True
        out = out and "workspace" in cont
        out = out and "version" in cont
        out = out and "license" in cont
        if not out:
            msg = "Invalid content of the application setup file ({}):\n{}"
            raise Exception(_Setup._setup_file, msg.format(cont))
        
    def reset_paths(self):
        """Reset workspace paths to package defaults.
        
        Define default paths for resource data and output folder.
        
        The default workspace is in user's Documents/stressfit folder. 
        If there is no user's Documents folder, use the current directory 
        as workspace.        
        """
        wks = self._data['workspace']
        docs = _Path.home().joinpath('Documents')
        if not docs.exists():
            wks['work'] = _Path().cwd()  
        else:
            wks['work'] = docs.joinpath('stressfit')
             
        wks['data'] = get_resource_path('data')
        wks['tables'] = get_resource_path('tables')
        wks['instruments'] = get_resource_path('instruments')
        wks['output']  = _Path('output') # relative to work


    def save(self):
        """Save application setup file."""
        try:
            wks = _workspace_to_str(self._data['workspace'])
            content = copy.deepcopy(self._data)
            content['workspace'].update(wks)
            txt = json.dumps(content, indent=4)
            with open(self._appsetup, 'w') as f:
                f.write(txt)
        except Exception as e:
            print('Cannot save application setup file: {}'.format(self._appsetup))
            print(e)

    def load(self):
        """Load application setup file."""
        try:
            with open(self._appsetup, 'r') as f:
                lines = f.readlines() 
            content = json.loads('\n'.join(lines))
            self.verify(content)
            wks = _str_to_workspace(content['workspace'])
            content['workspace'].update(wks)
            self._data.update(content)
        except Exception as e:
            print('Cannot load application setup file: {}'.format(self._appsetup))
            print(e)
            
    def info(self):
        fmt = '{}:\t{}'
        for key in self._data:
            print(fmt.format(key, self._data[key]))



class Workspace():
    """
    Class for managing user workspaces.
    
    A workspace is defined by the working root directory - `work`, and several
    other directories which help to structure input and output data:
    
    data:
        Input data (e.g. experimental data)
    tables:
        Lookup tables, material tables, etc.
    instruments:
        Instrument files (e.g. instrument configurations)
    output:
        Output data. All output files will be saved in this directory.
        
    If the other directories are given as relative paths, they are interpreted
    relative to the main `work` directory.
    
    The paths are internally stored as pathlib.Path objects.
    
    The workspace paths can be saved in the work directory in 
    the file `.stressfit/workspace.json` and loaded from it. 
    The file is searched for when switching to a new work directory 
    and loaded if found. 
    
    By default, `data`, `tables` and `instruments` point to the package 
    resources, `work` and `output` point to the user's Documents folder.
    
    Usage
    -----
    
    Use `change_workspace(root=newpath)` for switching to another root  
    directory.
    
    Use `set_paths(new_paths)` to change the workspace paths except root.
    
    """
    
    _wks_dir = '.stressfit' # directory with workspace info
    _wks_file = 'workspace.json' # file with workspace info
    def __init__(self):
        # set application default workspace
        self._logger = logger()
        self.reset_paths(create_tree=True, save=True)
    
    @property
    def wksdir(self):
        """Directory with workspace configuration data."""
        return self._paths['work']
    
    @property
    def keys(self):
        """Names of the workspace folder types.
        
        Example: 'work', 'data', 'tables', 'output', ...
        
        Used as the kind argument in various dataio functions.
        """
        return _wks_keys

    def _get_workspace_file(self):
        """Get workspace definition file.
        
        Create workspace definition folder if not present.
        """
        ep = self._paths['work']
        epath = _Path(ep).joinpath(Workspace._wks_dir)
        if not epath.exists():
            epath.mkdir(parents=True)
        elif epath.is_file():
            msg = 'Cannot create workspace folder. Remove file {} first.'
            raise Exception(msg.format(epath.as_posix()))
        # file with workspace info
        out = epath.joinpath(Workspace._wks_file)
        return out

    def _load(self):
        """Load workspace paths from local configuration file.
        
        Return
        ------
        dict
            Workspace paths. Return None if not loaded correctly.
        """
        wks = None
        try:
            file = self._get_workspace_file()
        except Exception as e:
            print(e)
            return None
        if file.exists():
            try:
                with open(file, 'r') as f:
                    lines = f.readlines() 
                res = json.loads('\n'.join(lines))
                if not 'workspace' in res:
                    msg = 'Cannot find workspace info in {}'
                    raise Exception(msg.format(file))
                lst = res['workspace']
                if not 'work' in lst:
                    lst['work'] = self._paths['work']
                wks = _str_to_workspace(lst)
            except Exception as e:
                print(e)
        return wks
        
    def reset_paths(self, create_tree=False, save=False):
        """Reset workspace paths to package defaults."""
        stp = _setup()
        self._paths = copy.deepcopy(stp.workspace)
        self.change_workspace(root=self._paths['work'], 
                              create_tree=create_tree,
                              save=save,
                              verbose=False)
        
    def change_workspace(self, root=None, create_tree=True, save=True, 
                         verbose=True):
        """Change workspace root directory.
        
        Process dependences:
            - Set the new workspace root directory
            - Try to load workspace.json from the workspace
            - If not loaded/not found, then make sure that relative paths
              exist in the new workspace. Save the workspace file. If not, 
              convert them to absolute paths using the previous workspace.
            - Create workspace tree if required

        Parameters
        ----------
        root : str
            The new workspace root directory.
            Must be a full path or None. If None, use current directory.
        create_tree : bool
            Create workspace tree
        save : bool
            Save workspace file if there is no one yet.
        verbose : bool
            Print info about new workspace.
            
        """
        err_msg = 'Cannot create a new workspace.'
        
        # set workspace path and verify
        if root is None:
            p = _Path.cwd()
        else:
            p = _Path(root)
        if not p.is_absolute():
            msg = 'Cannot use relative workspace directory: {}\n{}'
            raise Exception(msg.format(root,err_msg))        
        if (not p.exists()) or p.is_file():
            msg = 'Directory does not exist: {}\n{}'
            raise Exception(msg.format(root,err_msg))
        self._paths['work'] = p
        self._logger.clear()
        if verbose:
            self._logger.info('Workspace switched to {}'.format(self._paths['work']))
            #print('Workspace switched to {}'.format(self._paths['work']))        
        # try to load workspace file
        wks = self._load()
        
        # Cannot load workspace file?
        if wks is None:
            # make paths relative to the new workspace
            for key in _path_keys:
                p = self._paths[key]
                self._paths[key] = _mk_path_relative(self._paths['work'], p)
            # save workspace file
            if save:
                self.save()
            # NOTE: save also creates tree
            elif create_tree:
                _create_workspace_tree(self._paths) 
        # Workspace file loaded?
        else:
            if verbose:
                fn = self._get_workspace_file()
                self._logger.info('Loaded workspace setting from {}'.format(fn))
                # print('Loaded workspace setting from {}'.format(fn))
            self._paths.update(wks)
        
        # create subdirectories if needed
        if create_tree:
            _create_workspace_tree(self._paths) 
        

    def full_path(self, key, as_posix=False):
        """Return full path for given directory key.
        
        Parameters
        ----------
        key: str
            Directory key. One of the keys defined by `self.keys`.
            
        as_posix: bool
            If True, return as string, otherwise return as a pathlib.Path.
        """
        if key in self._paths:
            if self._paths[key].is_absolute():
                res = self._paths[key]
            else:
                res = self._paths['work'].joinpath(self._paths[key])
        else:
            msg = 'Unknown directory key ({}). Use one of {}.'
            raise Exception(msg,format(key, self.keys))
        if as_posix:
            return res.as_posix()
        else:
            return res
    
    def path(self, key=None, as_posix=False):
        """Return path name for given key. Default is workspace root."""
        if key in self._paths:
            p = self._paths[key]
        else:
            p = self._paths['work']
        if as_posix:
            return p.as_posix()
        else:
            return p
            
    def set_paths(self, **kwargs):
        """
        Set paths for data input and output folders.
        
        Paths with `None` value remain unchanged.
        
        Parameters
        ----------
        data: str
            Input data (measured strain, intensity, etc.).
        tables: str
            Other auxilliary input such as lookup tables or material data.
        instruments: str
            Instrument configuration files.
        output: str
            Output folder.
            
        If the given path is a subdirectory of the work path, then it is saved
        as relative.        
        """
        for key in kwargs:
            if key in _path_keys and kwargs[key] is not None:
                self._paths[key] = _mk_path_relative(self._paths['work'], 
                                                     _Path(kwargs[key]))
    
    def get_paths(self, absolute=False, keys=None, as_posix=True):
        """Return current posix path names as dict.
        
        Parameters
        ----------
        absolute: bool
            Return as absolute paths.
        keys: list
            If defined, return only listed directories. Otherwise
            return all workspace directories.
        as_posix : bool
            Return paths as strings
        
        Return
        ------
        dict
            Workspace paths as pathlib.Path objects or strings.    
            Keys are the directory types defined for the workspace.
        """
        env = {}
        if keys is None:
            keys = self.keys
        for key in keys:
            p = self._paths[key]
            if absolute:
                p = self.full_path(key, as_posix=as_posix)
            else:
                p = self._paths[key]
                if as_posix:
                    p = p.as_posix()
            env[key] = p
        return env
    
    def validate_paths(self):
        """Check that the workspace paths exist."""
        fmt = 'Path not found: {}'
        for key in self.keys:
            f = self.full_path(key, as_posix=False)
            if not f.exists():
                raise Exception(fmt.format(f.as_posix()))



    def save(self):
        """Save workspace paths in local configuration file."""
        file = self._get_workspace_file()
        keys = _path_keys
        try:
            _create_workspace_tree(self._paths)
            lst = _workspace_to_str(self._paths, keys=keys)
            out = {'workspace':lst}
            txt = json.dumps(out,indent=4)
            with open(file, 'w') as f:
                f.write(txt)
        except Exception as e:
            print(e)

    def info(self, absolute=False):
        """Print ifnormation on workspace paths."""
        fmt = '{}: {}'
        lst = _workspace_to_str(self._paths, absolute=absolute)
        for key in lst:
            print(fmt.format(key,lst[key]))
            
class Param(dict):
    """
    Class encapsulating a single instrument parameter.
    
    Implements variables as properties
    
    Parameters
    ----------
    key: str
        Parameter name.
    value: float
        Value.
    descr: str
        Brief description string.
    """ 
    
    def __init__(self, key:str, value:float, descr:str):
        dict.__init__(self, value=value, descr=descr)
        self.keystr = key
    
    @property
    def key(self):
        """Key string."""
        return self.keystr
    
    @property
    def value(self):  
        """Value of the parameter as `float`."""
        return self['value']
    
    @property
    def descr(self):  
        """Brief parameter description."""
        return self['descr']
    
    def json(self):
        """Return value and descr as JSON string."""
        h = {'value':self.value, 'descr':self.descr}
        ret = "{}".format(h)
        return ret
    
    def __str__(self):
        """Return string representation of Param."""
        return self.json()
         

class Table():
    """
    A simple class for 2D data with labeled columns and comments.
    
    Objects of this class are used to simplify I/O operations 
    with tabeled data stored as text files. 
    
    Attributes
    ----------
    file: str
        File name (filled when loaded).
    comments: str or list
        Comment lines.
    colHeaders: list
        Column headers.
    colHeaders: list
        Column headers.
    cells: numpy.array
        Table data as 2D numpy.array.
    size: tuple
        Table dimensions.
    
    """  
    
    def __init__(self):
        self.file = ''
        self.comments = []
        self.colHeaders = []
        self.rowHeaders = []
        self.cells = None
        self.size = (0, 0)
        self._hasRowHeaders = False
    
            
    def _headerstr(self):
        """Return header as a comment string."""
        s = ''
        nc = len(self.colHeaders)
        if (nc > 0):
            s += '# Header: '
            nc = len(self.colHeaders)
            for i in range(nc):
                s += '{}'.format(self.colHeaders[i])
                if (i==nc-1):
                    s += '\n'
                else:
                    s += '\t'
        return s
    
    def _rowstr(self, irow):
        """Return given row as a tab separated string."""
        s = ''
        lbl=''
        nc = len(self.rows[irow])
        if (len(self.rowHeaders)>irow):
            lbl = self.rowHeaders[irow]
            s += '{}\t'.format(lbl)
        for i in range(nc):   
            s += '{:g}'.format(self.rows[irow][i])
            if (i==nc-1):
                s += '\n'
            else:
                s += '\t'
        return s     
    
    def __str__(self):
        """Return table as a string with tab separated columns."""
        s = ''
        if (len(self.comments) > 0):
            for i in range(len(self.comments)):
                L = self.comments[i]
                if (not L.startswith('#')):
                    s += '# '
                s += '{}\n'.format(self.comments[i])
        s += self._headerstr()
        for i in range(len(self.rows)):
            s += self._rowstr(i)
        return s  


    def colIndex(self, label):
        """Return the index of a column with given label."""
        i = -1
        try:
            i = self.colHeaders.index(label)
            if (self._hasRowHeaders):
                i -= 1
        except:
            print('Warning: column label[{}] not found\n'.format(label))
        return i

    
    @property
    def rows(self):
        """Rows as a list of numpy.array objects."""
        return list(self.cells)

    @property
    def cols(self):
        """Columns as a list of numpy.array objects."""
        return list(self.cells.T)
       

### Path names handling

def set_path(**kwargs):
    """Shortcut for setting workspace directories.
    
    Calls `workspace().set_paths(**kwargs)`.
    """
    wks = workspace()
    wks.set_paths(**kwargs)

def get_resource_path(name):
    """
    Return absolute resource path for given resource folder.
    
    Allowed names are:
    data :
        for test input data 
    tables :
        for tabeled data supplied with the package
    instruments :
        for instrument configurations provided by the package
    conf :
        stressfit configuration files
    
    Return
    ------
    :class:`pathlib.Path`
    """
    p = _Path(__file__)
    if name in ['data','tables','instruments','conf']:
        return p.parent.joinpath('resources',name)
    else:
        raise Exception('Unknown resource name: {}.'.format(name))




def derive_filename(file, ext='', sfx=''):
    """
    Derive file name by changing extension and/or adding a suffix.
    
    - If `ext` is non-empty, remove any existing file extension and add
    the new one.   
    - If `sfx` is non-empty, add the give suffix to before the extension.
    
    Examples
    --------
    ``derive_filename('/my/long/path/file.txt', ext='dat', sfx='old')``
    
    returns '/my/long/path/file_old.dat' 
       
    Parameters
    ----------
     file: str
         Original file name (may be a full path)
     ext: str
         Extension to be added.
     sfx: str
         A suffix to be added at the end of filename (before extension).
         
    Returns
    -------
    :class:`pathlib.Path`
    """
    known_ext = ['.dat','.png','.txt','.csv','.xlsx','.doc','.docx']
    f = _Path(file)
    e = f.suffix
    path = f.parent
    base = f.stem   
    if e and e not in known_ext:
        base = ''.join([base, e])
        e=''
    if sfx:
        base = '_'.join([base, sfx])    
    if ext:
        fn = '.'.join([base, ext])
    else:
        if e:
            fn = ''.join([base, e])
        else:
            fn = base
    
    return path.joinpath(fn)


def get_input_file(filename, kind='data', path=None, **kwargs):
    """
    Convert filename to full path name.
    
    If path is not defined, use workspace path for given directory kind.
    Otherwise only convert the file name to a Path object.
    
    Parameters
    ----------
    filename: str
        File name (base name or full path).
    kind: str
        Which kind of directory to use:
            - `data` for input data
            - `tables` for table and other files
            - `instruments` for a file with instrument parameters
            - any other: current directory.
    path: str
        Optional search path for the input file. 
    
    Returns
    -------
    :obj:`pathlib.Path`
        Full path specification.
    """
    f = _Path(filename)
    if not f.is_absolute():
        wks = workspace()
        if path:
            p = _Path(path)
        elif kind in _wks_keys:
            p = wks.full_path(kind, as_posix=False)
        else:
            p = f.cwd()
        f = p.joinpath(f)
    return f


def get_output_file(filename, path=None, **kwargs):
    """
    Convert filename to full output path name.
    
    If path is not defined, use workspace path for output directory.
    Otherwise only convert the file name to a Path object.
    
    Parameters
    ----------
    filename: str
        File name (base name or full path).
    path: str
        Output path. If not defined and `filename` is relative, 
        then the default output path is used (see function :func:`.set_path`).
        
    Returns
    -------
    :class:`pathlib.Path`
        Full path specification as a :class:`pathlib.Path` object.
    """
    f = _Path(filename)
    if not f.is_absolute():
        if path:
            p = _Path(path)
        else:
            p = workspace().full_path('output', as_posix=False)
        f = p.joinpath(f)
    return f


### Decorators


def loadwrapper(func):
    """
    Decorate functions for loadig data.
    
    Create valid absolute path and print progress messages.
    
    On exception, print error message and re-raise the exception. 
    """
    @wraps(func)
    def func_wrapper(filename, **kwargs):
       fname = filename
       try:
           fname = get_input_file(filename, **kwargs)
           if fname.exists():
               #print('Trying to load {}.'.format(fname))
               res = func(fname, **kwargs)
               #print('Succeeded.')
           else:
               raise Exception('File {} does not exist.'.format(fname))
       except Exception as e:
           print('ERROR: could not load file {}.'.format(fname))
           print('Arguments: {}.'.format(kwargs))
           print('Check input path.')
           raise e
       return res
    return func_wrapper


def savewrapper(func):
    """
    Decorate functions for saving data.
    
    Create valid absolute path and print progress messages.
    
    On exception, print error message and re-raise the exception. 
    """
    @wraps(func)
    def func_wrapper(data, filename, **kwargs):
       fname = filename 
       try:
           fname = get_output_file(filename, **kwargs)
           func(data, fname, **kwargs)
           print('File {} saved.'.format(fname))
       except:
           print('Warning: could not save file {}.'.format(fname))
    return func_wrapper


### File read/write


def __std_header(comment="", source=None, labels=None):
    """Format standard header for file output.
    
    Parameters
    ----------
    comment: str or list
        Any comment lines.
    source: str
        Add 'Source: source' to the comments.
    labels: list
        Add a tab-separated list of labels.
        
    """
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = []
    out.append('Created: {}'.format(now))
    if source:
        out.append('Source: {}'.format(source))
    if comment:
        if (isinstance(comment,str)):
            out.append('{}'.format(comment))
        else:
            for i in range(len(comment)):
                out.append('{}'.format(comment[i]))   
    if labels:
        hdr = 'Header: ' + '\t'.join(labels)
        out.append(hdr)
    return out



def load_resource(filename, folder):
    """Read content of a resource file as a list of text lines.
    
    Parameters
    ----------
    filename: str
        base file name 
    folder: str
        resource sub-folder name
    """
    path = get_resource_path(folder)
    fn = path.joinpath(_Path(filename))
    f = open(fn, 'r')
    lines = f.readlines() 
    f.close()
    return lines


def load_config(filename, folder='conf', ext='json'):
    """Read a JSON configuration file file from resources.
    
    Parameters
    ----------
    filename: str
        base file name 
    folder: str
        resource sub-folder name
    
    Returns
    -------
    dict
        The file content converted to a dictionary.
    """
    if not filename.endswith('.'+ext):
        filename += '.{}'.format(ext)
    conf = load_resource(filename, folder)
    try:
        res = json.loads('\n'.join(conf))
    except Exception as e:
        print(e)
    return res

@loadwrapper
def load_text(filename, **kwargs):
    """Load a text file as is.
                        
    Parameters
    ----------
    filename: str
        Input file name (base name or absolute path).
    **kwargs:
        Keyword parameters passed to :func:`.get_input_file`. Use ``kind`` to
        define default input path and ``path`` to override the default input 
        path. 

    Use :func:`~.set_path` to define default input/output folders.
    
    Returns
    -------
    :obj:`list` of :obj:`str`
        Content of the file as a list of lines.
    """
    f = open(filename, 'r')
    lines = f.readlines() 
    f.close()
    return lines


@savewrapper
def save_text(text, filename, **kwargs):
    """Save a text in a text file.
                        
    Parameters
    ----------
    text: :obj:`str` or :obj:`list` of :obj:`str`
        Text to save, either a string or a list of strings.
    filename: str
        Output file name (base name or absolute path).
    **kwargs:
        Keyword parameters passed to :func:`.get_output_file`. 
        Use ``path`` to override the default input path.
    
    Returns
    -------
    :obj:`list` of :obj:`str`
        Content of the file as a list of lines.
    """
    f = open(filename, 'w')
    if isinstance(text, list):
        out = '\n'.join(text)
        f.write(out+'\n')
    else:
        f.write(text+'\n')
    f.close()


@loadwrapper
def load_data(filename, rows=[0, -1], verbose=True, **kwargs):
    """
    Load a numeric table from a text file using :meth:`numpy.loadtxt`.
   
    Parameters
    ----------
    filename : str
        File name.
    rows : list(2)
        Range of rows (excluding # comments) to evaluate. 
        Use [0, -1] for full range.
    verbose : bool
        If True, print info on loaded file.
    **kwargs
        Keyword arguments passed to :func:`numpy.loadtxt`
    
    Returns
    -------
    :obj:`numpy.ndarray`
        Specified range of rows from the input file.
    """
    keys_to_skip = ['path', 'kind', 'rows']
    for k in keys_to_skip:
        if k in kwargs:
            kwargs.pop(k)
    d = np.loadtxt(filename, **kwargs)
    if verbose:
        print('File loaded: {}'.format(filename))
    if (rows[1]>=0):
        return d[rows[0]:rows[1]+1,:]
    elif rows:
        return d[rows[0]:,:]
    else:
        return d[:,:]


def load_table(filename, **kwargs):
    """
    Load a table with a single line header from a text file.
    
    The header signature is 'Header:' followed by the list of column labels.
    Except of the header, skip comment lines (#) and empty rows. 
                                              
    The table body must contain the same number of data columns as is the 
    number of header labels. Otherwise an exception is raised.
    
    The column separator is any whitespace character.
                                              
    Parameters
    ----------
    filename: :obj:`pathlib.Path`
        File name.
    **kwargs:
        Keyword parameters passed to :func:`.load_text`
            - Use ``kind`` to choose the default input path 
              (see :func:`~.set_path`).
            - Use ``path`` to override the default input path.
    
    Returns
    -------
    :class:`.Table`
    
    """
    ERR_STRUCT = 'Wrong table structure on line {}\n'
    ERR_COLNUM = 'Unequal number of columns, given {:d} expected {:d}\n'
    ERR_ROW_LABELS = 'Nonempty row labels are required at all rows or not at all\n'
    ERR_CELL_NUMBER = 'All rows must have the same number of cells\n'
    
    iline = 0
    ncols = 0
    nrows = 0
    def getHdr(line):
        keys = []
        key = 'Header:'
        i = line.find(key)
        if (i>0):
            j = i + len(key)
            lbls = line[j:].split()
            nc = len(lbls)
            for i in range(nc):
                keys.append(lbls[i])
        return keys
    
    def getRow(line):
        row = []
        items = line.split()
        nc = len(items)
        for i in range(nc):
            s = items[i].strip()
            t = -1
            try:
                x = int(s)
                row.append(x)
                t = 1
            except:
                try:
                    x = float(s)
                    row.append(x)
                    t = 2
                except:
                    row.append(s)
                    t = 0
            # only the 1st item can be string = row label
            if(i>0 and t<1):
                msg = 'Wrong table syntax at line[{:d}]\n'.format(iline)
                msg += 'All columns except of the row header must be numeric\n'
                raise Exception(msg)
        return row

    hdr = []
    table = None
    hasRowHeaders = False
    lines = load_text(filename,  **kwargs)
    table = Table()
    rows = []
    if isinstance(filename, _Path):
        table.file = filename.as_posix()
    else:
        table.file = filename
    for line in lines:
        iline +=1
        L = line.strip()
        # header or comment
        if (L.startswith('#')):
            # try header
            hdr = getHdr(L)
            if (ncols<1 and len(hdr)>0):
                table.colHeaders = hdr
                # number of columns is defined by the headder
                ncols = len(hdr)
            # anything else starting with # is a comment
            else:
                s = L.split("#")
                table.comments.append(s[1].strip())
            
        # table row
        elif len(L)>0:
            row = getRow(L)
            nc = len(row)
            # ignore empty rows
            if nc>0:
                # check number of cells in the row
                if (ncols<1):
                    # number of columns is defined by the 1st row
                    ncols = nc
                elif (ncols!=nc):
                    msg = ERR_STRUCT + ERR_COLNUM
                    raise Exception(msg.format(iline, nc, ncols))
                # pop a label if any                   
                if (isinstance(row[0],str)):
                    lbl = row.pop(0)
                    if (nrows==0):
                        hasRowHeaders = True
                    elif (not hasRowHeaders):
                        msg = ERR_STRUCT + ERR_ROW_LABELS
                        raise Exception(msg.format(iline))
                    table.rowHeaders.append(lbl)
                elif (hasRowHeaders):
                    # if the table includes row geaders, they must start every row
                    msg = ERR_STRUCT + ERR_ROW_LABELS
                    raise Exception(msg.format(iline))
                if (nrows>0):
                    # check that each row has the same number of elements
                    nr = len(rows[0])
                    if (len(row) != nr):
                        msg = ERR_STRUCT + ERR_CELL_NUMBER
                        raise Exception(msg.format(iline))
                rows.append(row)
                nrows += 1
    # consolidate table
    table._hasRowHeaders = hasRowHeaders
    if (nrows>0):
        ncols = len(rows[0])
        table.size = (nrows, ncols)
        table.cells = np.array(rows)
    else:
       print('Warning: the table {} is empty'.format(table.file))
    return table


def load_params(filename, **kwargs):
    """Load parameters from a text file.
    
    Parameters in the file must have the format ``key=value # description``.
    Empty lines and lines starting with # are ignored.
   
    The parameters are returned as a :obj:`dict` of :class:`~.Param`. 
    
    Parameters
    ----------
    filename: str
        File name.
    **kwargs:
        Keyword parameters passed to :func:`.load_text`
            - Use ``kind`` to choose the default input path 
              (see :func:`~.set_path`).
            - Use ``path`` to override the default input path. 

    Returns
    -------
    :obj:`dict` of :class:`~.Param`    

    """
    fmtw = "WARNING: Duplicite key {} found. Using the last definition."
    h = {}
    lines = load_text(filename,  **kwargs)
    for ln in lines:
        L = ln.strip()
        if (len(L)>0) and not L.startswith('#'):
            row = L.split('=')
            if (len(row)>1):
                key = row[0].strip()
                if key in h:
                    print(fmtw.format(key))
                rest = row[1].split('#')
                value = rest[0].strip()
                if len(rest)>1:
                    descr = rest[1].strip()
                else:
                    descr = ''
                fvalue = float(value)
                p = Param(key, fvalue, descr)
                h[key] = p
            else:
                msg = 'Invalid format, expected [key=value # description].'
                raise Exception('{}\n'.format(L)+msg)
    return h            


@savewrapper
def save_data(data, filename, header=[], comment="", source="", **kwargs):
    """
    Save an array as a text table using :func:`numpy.savetxt`.
    
    Parameters
    ----------
    data: array or list
        Data to be saved.
    filename: str
        Filename for output. 
        If path is not included, the default output path is used. 
        Use :func:`~.set_path` to define the output folder.
    header: list
        List of column labels.
    comment: str or list of str
        Header comments.
    source: str
        Optionally, add id of the calling script to the coments.
        Use "``source=__file__``" to add the name of the calling script)
    **kwargs:
        Keyword parameters passed to :func:`~.get_output_file`. 
        Use ``path`` to override the default output path.
        
    """
    comstr = ""
    hdr = __std_header(comment=comment, source=source, labels=header)
    for h in hdr:
        comstr += '{}\n'.format(h)
    np.savetxt(filename, np.array(data), delimiter='\t', header=comstr, 
               fmt="%g")


def save_table(data:Table, filename, source="", **kwargs):
    """
    Save Table object as a text file.
    
    Parameters
    ----------
    data: :class:`~.Table`
        Table to be saved.
    filename: str
        Filename for output. For relative paths, the default output path 
        is used. Use :func:`~.set_path` to define the output folder.
    source: str
        Source name passed to :func:`.save_data`.
    **kwargs:
        Keyword parameters passed to :func:`~.get_output_file`. 
        Use ``path`` to override the default output path.        
        
    """
    hdr = __std_header(comment=data.comments, 
                       source=source, 
                       labels=data.colHeaders)
    out = []
    nr = data.size[0]
    nc = data.size[1]
    for h in hdr:
        out.append('# {}'.format(h))
    if len(data.rowHeaders)==nr:
        fmt = '{}' + nc*'\t{:.10g}'
        for i in range(data.size[0]):
            ln = fmt.format(data.rowHeaders[i], *data.cells[i,:])
            out.append(ln)
    else:
        fmt = '{:g}' + (nc-1)*'\t{:.10g}'
        for i in range(nr):
            ln = fmt.format(*data.cells[i,:])
            out.append(ln)
    save_text(out, filename, **kwargs)


def save_params(data, filename, comment="", source="", **kwargs):
    """
    Save parameters as a text file.
    
    Parameters
    ----------
    data: :obj:`dict` of :class:`~.Param`
        Parameters represented by the  :class:`.Param` objects.
    filename: str
        Filename for output. For relative paths, the default output path 
        is used. Use :func:`~.set_path` to define the output folder.
    comment: :obj:`str` or :obj:`list` of :obj:`str`
        Header comments.
    source: str
        Optionally, add id of the calling script to the comments.
        Use "``source=__file__``" to add the name of the calling script)
    **kwargs:
        Keyword parameters passed to :func:`~.save_text`. 
        Use ``path`` to override the default output path.
        
    """
    out = []
    hdr = __std_header(comment=comment, source=source)
    for h in hdr:
        out.append('# {}'.format(h))
    fmt1 = '{}={:g} # {}'
    fmt2 = '{}={:g}'
    for key in data:
         p:Param = data[key]
         if p.descr:
             s = fmt1.format(key,p.value, p.descr)
         else:
             s = fmt2.format(key,p.value)
         out.append(s)
    save_text(out, filename, **kwargs)




def test():
    """Test unit."""
    # wks = workspace()
    # wks.info()
    stp = _setup()
    stp.info()
    print()
    wks = workspace()
    wks.info()
    d = load_data('eps_B_axi.dat', kind='data')
    assert len(d)>10


#log = logging.getLogger()
#try: 
#    a=1/0
#except Exception as e: 
#    log.exception('aaa')

    