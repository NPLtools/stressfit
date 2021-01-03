# -*- coding: utf-8 -*-
# Created on Mon Jun  4 10:50:25 2018
# @author: User
"""
Input/output methods for data used by StressFit.
- Based on the `pathlib` package.
- Implements a simple Table class for 2D data arrays with headers and comments.
- Handles default repository paths


"""
import numpy as np
import datetime
from pathlib import Path


def get_resource_path(name):
    """
    Return absolute resource path for given resource folder.
    Allowed names are:
        - `data` for test input data 
        - `tables` for tabeled data supplied with the package
        - `instruments` for instrument configurations provided by the package
    """
    
    p = Path(__file__)
    return p.parent.joinpath('resources',name)


# Define default paths for resource data and output folder.
# These values can be overriden by set_path().
__inpath = get_resource_path('data')
__tables = get_resource_path('tables')
__instruments = get_resource_path('instruments')
__outpath = Path().cwd()


def set_path(inpath=None, tables=None, instruments=None, outpath=None):
    """
    Set paths for data input and outpout folders. 
    By default, the input paths are the package resource directories, 
    the output path is the current directory. 
    Paths with `None` value remain unchanged.
    
    Parameters
    ----------
    inpath: str
        Input / experimental data.
    tables: str
        Other auxilliary input such as lookup tables or material data.
    instruments: str
        Instrument configuration files
    outpath: str
        Output folder
    
    """    
    
    global __inpath, __tables, __instruments, __outpath
    if (inpath):
        __inpath = Path(inpath)
    if (tables):
        __tables = Path(tables)
    if (instruments):
        __instruments = Path(instruments)
    if (outpath):
        __outpath = Path(outpath)
        

### Classes


class Table():
    """
    A simple class for 2D data with labeled columns and comments.
    
    Objects of this class are used to simplify I/O operations 
    with tabeled data stored as text files. 
    
    """  
    
    def __init__(self):
        self.file = ''
        self.comments = []
        self.colHeaders = []
        self.rowHeaders = []
        self.cells = None
        self.rows = []
        self.columns = []
        self.size = (0, 0)
        self.hasRowHeaders = False
    
    def colIndex(self, label):
        """Return the index of a column with given label."""
        i = -1
        try:
            i = self.colHeaders.index(label)
            if (self.hasRowHeaders):
                i -= 1
        except:
            print('Warning: column label[{}] not found\n'.format(label))
        return i
            
    def headerStr(self):
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
    
    def rowStr(self, irow):
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
    
    def toText(self):
        """Return table as a string with tab separated columns."""
        s = ''
        if (len(self.comments) > 0):
            for i in range(len(self.comments)):
                L = self.comments[i]
                if (not L.startswith('#')):
                    s += '# '
                s += '{}\n'.format(self.comments[i])
        s += self.headerStr()
        for i in range(len(self.rows)):
            s += self.rowStr(i)
        return s   


### Path names handling

def derive_filename(file, ext='', sfx=''):
    """Derive file name by changing extension and/or adding a suffix. 
    
    - If `ext` is non-empty, remove any existing file extension and add
    the new one.   
    - If `sfx` is non-empty, add the give suffix to before the extension.
    
    Example
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
         
    
    """
    
    f = Path(file)
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


def get_input_file(fname, kind='input'):
    """ Convert fname to full path name. 
    If path is not included, prepend the default input path (see set_path).
    Otherwise only convert the file name to a Path object.
    
    Parameters
    ----------
    
    fname: str
        File name (base or full path).
    kind: str
        Search path specification:
            - `input` for input data
            - `table` for table and other files
            - `instrument` for a file with isntrument parameters
            - any other: current directory.
    
    
    Returns
    -------
    
    Path
        Full path specification as a pathlib.Path object.
    """
    
    f = Path(fname)
    if not f.is_absolute():
        if (kind == 'input'):
            p = __inpath
        elif (kind == 'table'):
            p = __tables
        elif (kind == 'instrument'):
            p = __instruments
        else:
            p = f.cwd()
        f = f.joinpath(p,f)
    return f


def get_output_file(fname):
    """ Convert fname to full path name.
    If path is not included, prepends the default output path (see set_path).
    """
    f = Path(fname)
    if not f.is_absolute():
        f = f.joinpath(__outpath,f)
    return f

### Decorators


def loadwrapper(func):
    """ Decorator applied on functions loadig input data. On exception, it 
    prints error message and re-raises the exception. 
    """
    
    def func_wrapper(name, kind, **kwargs):
       try:
           fname = get_input_file(name, kind=kind)
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
    """ Decorator applied on functions saving auxilliary data.
    Calls func(name) and prints a message if successful.
    On exception prints a warning message. 
    """
    
    def func_wrapper(data, name, **kwargs):
       fname = name 
       try:
           fname = get_output_file(name)
           func(data, fname, **kwargs)
           print('File {} saved.'.format(fname))
       except:
           print('Warning: could not save file {}.'.format(fname))
    return func_wrapper



### File read/write

@loadwrapper
def load_data(filename, kind='input', rmin=0, rmax=-1, **kwargs):
    """Read a numeric table from given text file using numpy.loadtxt().
    Return specified range of rows.
    
    
    Parameters
    ----------
    
    filename: str
        File name.
    rmin, rmax: int
        Range of rows (excluding # comments) to evaluate. Use rmax=-1 for 
        full range.
    kind: str
        A string passed to get_input_file as the `kind` parameter value.
        It determines the search path for the data. The allowed values are:
            - `input` for input data
            - `table` for table and other files
            - `instrument` for a file with isntrument parameters
            - empty string: current directory.  
    
    kwargs
        Keyword arguments passed to numpy.loadtxt
    
    Returns
    -------
    
    ndarray
        Specified range of rows from the input file.
    
    Note
    ----
    
    Use :func:`~stressfit.dataio.set_path` to define default input/output folders.
    """
    d = np.loadtxt(filename, **kwargs)
    if (rmax>=0):
        return d[rmin:rmax+1,:]
    else:
        return d[rmin:,:]

@loadwrapper
def read_dict(filename, kind='input'):
    """Reads a simple dictionary in the format `name=value`.
    Skip comment lines (#). 
                        
    Parameters
    ----------
    
    filename: str
        File name.
    kind: str
        A string passed to get_input_file as the `kind` parameter value.
        It determines the search path for the data. The allowed values are:
            - `input` for input data
            - `table` for table and other files
            - `instrument` for a file with isntrument parameters
            - empty string: current directory.            
    
    Use :func:`~stressfit.dataio.set_path` to define default input/output folders.
    
    Returns
    -------
    dict
        A hash map with keys and values found in the file.
    
    """
    
    h = {}
    fname = get_input_file(filename, kind=kind)
    if (fname.exists()):
        f = open(filename, 'r')
        lines = f.readlines() 
        f.close()
        for line in lines:
            L = line.strip()
            if (not L.startswith('#')):
                row = L.split('=')
                if (len(row)>1):
                    key = row[0].strip()
                    value = row[1].split('#')[0].strip()
                    try:
                        h[key] = float(value)
                    except:
                        print("ERROR: Can''t read parameter {} = {}".format(key,value))
    return h

@loadwrapper
def read_table(filename, kind='table'):
    """ Read a table with a single line header.
    The header signature is 'Header:' followed by the list of column labels.
    Except of the header, skip comment lines (#) and empty rows. 
                                              
    The table body must contain the same number of data columns as is the 
    number of header labels. Otherwise an exception is raised.
    
    The column separator is any whitespace character.
                                              
    Parameters
    ----------
    
    filename: str
        File name.
    kind: str
        A string passed to get_input_file as the `kind` parameter value.
        It determines the search path for the data. The allowed values are:
            - `input` for input data
            - `table` for table and other files
            - `instrument` for a file with isntrument parameters
            - empty string: current directory.            
    
    Use :func:`~stressfit.dataio.set_path` to define default input/output folders.
    
    Returns
    -------
    stressfit.Table
    
    """
    ERR_STRUCT = 'Wrong table structure on line {}\n'
    ERR_COLNUM = 'Unequal number of columns, expected [{:d}]\n'
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
    fname = filename
#    fname = get_input_file(filename, kind=kind)
    if (fname.exists()):
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
        table = Table()
        table.file = fname.as_posix()
        for line in lines:
            iline +=1
            L = line.strip()
            # header or comment
            if (L.startswith('#')):
                # try header
                hdr = getHdr(line)
                if (ncols<1 and len(hdr)>0):
                    table.colHeaders = hdr
                    # number of columns is defined by the headder
                    ncols = len(hdr)
                # anything else starting with # is a comment
                else:
                    table.comments.append(L)
                
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
                        raise Exception(msg.format(iline,ncols))
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
                        msg = ERR_STRUCT + ERR_ROW_LABELS
                        raise Exception(msg.format(iline))
                    if (nrows>0):
                        nr = len(table.rows[0])
                        if (len(row) != nr):
                            msg = ERR_STRUCT + ERR_CELL_NUMBER
                            raise Exception(msg.format(iline))
                    table.rows.append(row)
                    nrows += 1
        # consolidate table
        table.hasRowHeaders = hasRowHeaders
        if (nrows>0):
            ncols = len(table.rows[0])
            table.size = (nrows, ncols)
            table.cells = np.array(table.rows)
            table.columns = []
            for i in range(ncols):
                table.columns.append(table.cells[:,i])
        else:
           print('Warning: the table {} is empty'.format(table.file))
    return table


@savewrapper
def save_data(data, filename, header=[], comment="", source=""):
    """Save an array as a text table.
    
    Parameters
    ----------
    
    data: array or list
        Data to be saved
    filename: str
        Filename for output. If path is not included, the default output 
        path is used. 
        Use :func:`~stressfit.dataio.set_path` to define the output folder.
    header: list
        List of column labels
    comment: str or list of str
        Header comments
    source: str
        optional id of the calling script (use source=__file__)
        
    """
    
    
    now = datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
    comstr = 'Created: {}\n'.format(now)
    if (source):
        comstr += 'Source: {}\n'.format(source)
    if (isinstance(comment,str)):
        comstr += '{}\n'.format(comment)
    else:
        for i in range(len(comment)):
            comstr += '{}\n'.format(comment[i])
    hdr = "Header:"
    for i in range(len(header)):
        hdr += '\t{}'.format(header[i])
    np.savetxt(filename, np.array(data), delimiter='\t', 
               header=comstr+hdr, fmt="%g")


def save_table(data, filename, source=""):
    """Save Table object as a text file.
    
    Parameters
    ----------
    
    data: Table
        Data to be saved
    filename: str
        Filename for output. If path is not included, the default output 
        path is used. 
        Use :func:`~stressfit.dataio.set_path` to define the output folder.
        
    """
    save_data(data.cells, filename, header=data.colHeaders, 
              comment=data.comments, source=source)

    