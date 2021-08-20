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
from pathlib import Path as _Path
from functools import wraps


### Classes


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
        - `data` for test input data 
        - `tables` for tabeled data supplied with the package
        - `instruments` for instrument configurations provided by the package
    
    Return
    ------
    :class:`pathlib.Path`
    """
    p = _Path(__file__)
    if name in ['data','tables','instruments']:
        return p.parent.joinpath('resources',name)
    else:
        raise Exception('Unknown resource name: {}.'.format(name))


def set_path(data=None, tables=None, instruments=None, output=None):
    """
    Set paths for data input and outpout folders.
    
    By default, the input paths are the package resource directories, 
    the output path is the current directory. 
    Paths with `None` value remain unchanged.
    
    Parameters
    ----------
    data: str
        Input of strain and intensity data.
    tables: str
        Other auxilliary input such as lookup tables or material data.
    instruments: str
        Instrument configuration files.
    output: str
        Output folder.
    
    """
    global __inpath, __tables, __instruments, __outpath
    if (data):
        __inpath = _Path(data)
    if (tables):
        __tables = _Path(tables)
    if (instruments):
        __instruments = _Path(instruments)
    if (output):
        __outpath = _Path(output)


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
    Convert fname to full path name.
    
    If path is not included, prepend the default input path (see set_path).
    Otherwise only convert the file name to a Path object.
    
    Parameters
    ----------
    filename: str
        File name (base name or full path).
    kind: str
        Search path specification:
            - `data` for input data
            - `tables` for table and other files
            - `instruments` for a file with instrument parameters
            - any other: current directory.
    path: str
        Optional search path for the input file. 
        If defined and if ``filename`` is
        a relative path, it overrides the default search path (see 
        the parameter ``kind`` and function :func:`.set_path`).
    
    Returns
    -------
    :class:`pathlib.Path`
        Full path specification as a :class:`pathlib.Path` object.
    """
    f = _Path(filename)
    if not f.is_absolute():
        if path:
            p = path
        elif (kind == 'data'):
            p = __inpath
        elif (kind == 'tables'):
            p = __tables
        elif (kind == 'instruments'):
            p = __instruments
        else:
            p = f.cwd()
        f = f.joinpath(p,f)
    return f


def get_output_file(filename, path=None, **kwargs):
    """
    Convert fname to full path name.
    
    Parameters
    ----------
    fname: str
        File name (base name or full path).
    path: str
        Output path. If not defined and ``filename`` is relative, 
        then the default output path is used (see function :func:`.set_path`).
        
    Returns
    -------
    :class:`pathlib.Path`
        Full path specification as a :class:`pathlib.Path` object.
    """
    f = _Path(filename)
    if not f.is_absolute():
        if path:
            f = f.joinpath(path,f)
        else:
            f = f.joinpath(__outpath,f)
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


# Define default paths for resource data and output folder.
# These values can be overriden by set_path().
__inpath = get_resource_path('data')
__tables = get_resource_path('tables')
__instruments = get_resource_path('instruments')
__outpath = _Path().cwd()


