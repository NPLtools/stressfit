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
from pathlib import Path as _Path
from functools import wraps

__work = None



### Classes

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
    
    The paths are internally saved as pathlib.Path objects.
    
    The workspace paths can be saved in the work directory in 
    the file `.stressfit` and loaded from it. The file is searched for when
    switching to a new work directory and loaded if found. 
    
    By default, `data`, `tables` and `instruments` point to the package 
    resources, `work` and `output` point to the user's Documents folder. 
    
    
    """
    cfg_name = '.stressfit'
    types = ['work', 'data', 'tables', 'instruments', 'output']
    def __init__(self):
        self._path_keys = Workspace.types
        self._paths = {} 
        self.reset_paths()

    def reset_paths(self):
        """Reset paths to package defaults.
        
        Define default paths for resource data and output folder.
        These values can be overriden by set_path().
        """
        self._paths['data'] = get_resource_path('data')
        self._paths['tables'] = get_resource_path('tables')
        self._paths['instruments'] = get_resource_path('instruments')
        self._paths['output']  = _Path.home().joinpath('Documents')
        self._paths['work'] = self._paths['output']
        if not self._paths['output'].exists():
            self._paths['output'] = _Path().cwd()
            self._paths['work'] = self._paths['output']
  
    def keys(self):
        return self._path_keys
  
    def full_path(self, key, as_posix=False):
        """Return full path for given directory key.
        
        Parameters
        ----------
        key: str
            Directory key. One of the keys defined by `self._path_keys`.
            
        as_posix: bool
            If True, return as string, otherwise return as a pathlib.Path.
        """
        if key in self._path_keys:
            if key=='work':
                res = self._paths['work']
            elif self._paths[key].is_absolute():
                res = self._paths[key]
            else:
                res = self._paths['work'].joinpath(self._paths[key])
        else:
            msg = 'Unknown directory key ({}). Use one of {}.'
            raise Exception(msg,format(key, self._path_keys))
        if as_posix:
            return res.as_posix()
        else:
            return res
    

    def change_workspace(self, work_path):
        """
        Change workspace root directory.
        
        Process dependences:
            
            - Try to load the local .stressfit configuration file
            - If not loded/not found, then make sure that relative paths
              exist in the new workspace. 
            - If not, convert them to absoluet paths.

        Parameters
        ----------
        work_path : str
            Must be a full path. If not, :meth:`pathlib.Path.absolute` is 
            tried to derive absolute path.
        
        Return
        ------
        bool
            True if new workspace configuration was loaded. 
        """
        out = False
        old_work = self._paths['work']
        p = _Path(work_path)
        if p.is_absolute():
            self._paths['work'] = p
        else:
            try:
                pp = p.absolute()
                self._paths['work'] = pp                
            except Exception:
                msg = 'Cannot set relative work path: {}'
                print(msg.format(work_path))
        # try to load configuration
        out = self.load()
        # if not loaded, check relative paths:
        # - keep it relative if it exists in the new wrokspace
        # - otherwise convert to absolute
        if not out:
            for key in self._paths.keys():
                if key == 'work':
                    pass
                else:
                    p = self._paths[key]
                    if not p.is_absolute():
                        pp = self._paths['work'].joinpath(p)
                        if not pp.exists():
                            self._paths[key] = old_work.joinpath(p)
            self.set_paths()
        return out
        
        
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
            if key == 'work':
                pass
            elif key in self._path_keys and kwargs[key] is not None:
                p = _Path(kwargs[key])
                # make absolute paths relative to work if possible
                if p.is_absolute():
                    try:
                        p = p.relative_to(self._paths['work'])
                    except:
                        pass    
                self._paths[key] = p
    
    def get_paths(self, absolute=False, keys=None):
        """Return current posix path names as dict.
        
        Parameters
        ----------
        absolute: bool
            Return as absolute paths.
        keys: list
            If defined, return only listed directories. Otherwise
            return all workspace directories.
        
        Return
        ------
        dict
            Keys are the directory types defined for the workspace.
            The values are path names returned by :func:`pahlib.Path.as_posix`.
        """
        env = {}
        if keys is None:
            keys = self._path_keys
        for key in keys:
            p = self._paths[key]
            if absolute:
                p = self.full_path(key, as_posix=True)
            else:
                p = self._paths[key].as_posix()
            env[key] = p
        return env
    
    def validate_paths(self):
        """Check that the workspace paths exist."""
        fmt = 'Path not found: {}'
        for key in self._path_keys:
            f = self.full_path(key, as_posix=False)
            if not f.exists():
                raise Exception(fmt.format(f.as_posix()))

    def load(self):
        """Load workspace path names from local workspace configuration file.
        
        If successful, update workspace setting.
        
        Return
        ------
        bool
            True if successfully loaded and workspace updated.
        """
        out = False
        file = self._paths['work'].joinpath(Workspace.cfg_name)
        if file.exists():
            try:
                f = open(file, 'r')
                lines = f.readlines() 
                f.close()
                res = json.loads('\n'.join(lines))
                if 'workspace' in res:
                    w = res['workspace']
                    self.set_paths(**w)
                    out = True
            except Exception as e:
                print(e)
        return out

    def save(self):
        """Save workspace path names in local workspace configuration file."""
        if self._paths['work'].exists():
            file = self._paths['work'].joinpath(Workspace.cfg_name)
            keys = self._path_keys.copy()
            try:
                keys.remove('work')
                lst = self.get_paths(absolute=False, keys=keys)
                out = {'workspace':lst}
                txt = json.dumps(out,indent=4)
                f = open(file, 'w')
                f.write(txt)
                f.close()
            except Exception as e:
                print(e)

    def print_info(self, absolute=False):
        """Print ifnormation on workspace paths."""
        fmt = '{}: {}'
        lst = self.get_paths(absolute=absolute)
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
    f = _Path(file)
    e = f.suffix
    path = f.parent
    base = f.stem       
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
        if path:
            p = _Path(path)
        elif kind in __work.keys():
            p = __work.full_path(kind, as_posix=False)
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
            p = __work.full_path('output', as_posix=False)
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
               print('Trying to load {}.'.format(fname))
               res = func(fname, **kwargs)
               print('Succeeded.')
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
def load_data(filename, rows=[0, -1], **kwargs):
    """
    Load a numeric table from a text file using :meth:`numpy.loadtxt`.
   
    Parameters
    ----------
    filename: str
        File name.
    rows: list(2)
        Range of rows (excluding # comments) to evaluate. 
        Use [0, -1] for full range.
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
        Optionally, add id of the calling script to the coments.
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


def test():
    w = workspace()
   #  w.print_info(absolute=True)
    d = load_data('eps_B_axi.dat', kind='data')
    assert len(d)>10

    
    
