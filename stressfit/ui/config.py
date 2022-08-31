# -*- coding: utf-8 -*-
"""Configuration and input data for user interface.

Created on Mon Apr 25 15:47:30 2022
@author: Jan Saroun, saroun@ujf.cas.cz
"""

# TODO connect with a messaging tool instead of print()
# TODO define strain distribution list, linked to input data
# TODO define command for plotting & fitting strain  model + data)

import stressfit.dataio as dataio
import stressfit.commands as comm
from stressfit.geometry import Geometry
import copy
import re
import json

__config = None

def uiconfig():
    """Access UI_config object.
    
    UI_config is created when this function is called first time and persists
    as a singleton in given kernel scope.
    """
    global __config
    if __config is None:
        __config = UI_config()
    return __config

class ScanData():
    """Encapsulate input (experimental) data.
    
    Use to load strain and intensity data and assign geometry and sampling.
    
    Properties
    ----------
    scan: dict
        Returns scan data with assigned geometry and sampling data, 
        compatible with :func:`stressfit.commands.set_scan`.
    
    """
    
    def __init__(self):
        self._scan = {}
       
    @property
    def scan(self):
        """Scan data with assigned geometry and sampling data.
        
        Compatible with :func:`stressfit.commands.set_scan`.
        """
        return self._scan
    
    def load(self, fstrain, fint=None, path="", verbose=False):
        """Load scan experimental data.

        Parameters
        ----------
        fstrain : str
            File name for strain scan data (pos, eps, error).
        fint : str, optional
            File name for intensity scan data (pos, intensity, error).
        path : str, optional
            Input path. If not defined, file is searched in the workspace 'data' folder.
        verbose : bool, optional
            Print information on file load.

        """
        self._scan['eps'] = dataio.load_data(fstrain, path=path, verbose=verbose)
        self._scan['epsfile'] = fstrain
        if fint:
            try:
                self._scan['int'] = dataio.load_data(fint, path=path, verbose=verbose)
                self._scan['intfile'] = fint
            except:
                self._scan['int'] = None
                self._scan['intfile'] = ''
        else:
            self._scan['int'] = None
            self._scan['intfile'] = ''

    def assign_geometry(self, geom):
        """Assign geometry to the scan data.

        Parameters
        ----------
        geom : dict
            Geometry vectors - see :class:`stressfit.geometry.Geometry`.

        """
        scan_geometry = {key: geom[key] for key in Geometry.input_keys}
        self._scan.update(**scan_geometry)
    
    def assign_sampling(self, sampling):
        """Assign sampling data to the scan.
        
        Parameters
        ----------
        sampling : obj
            Sampling data object as returned 
            by :meth:`stressfit.commands.load_sampling`.

        """
        self._scan['sampling'] = sampling
        
    

class UI_data():
    """Container for user input data.
    
    Initial content is loaded from resources (config/config_ui.json).   
    
    Parameters
    ----------
    config : dict
        Configuration parameters to get initial data from.
        It must include the keys `udata` and `uinput` with corresponding
        user data and user input.
    
    Data names
    ----------
    
    shape: dict
        Sample shape parameters       
    geometry: list
        List of geometry parameters
    sampling: list
        List of sampling data sets
    attenuation: dict
        Attenuation data
    data: list
        List of experimental data sets 
        (intensity, strain and associated geometry and sampling)
    imodel: list
        List of intensity models
    emodel: list
        List of strain models
    
    
    Usage
    -----
    load(name, item=...)
        Load required file input.
    reload_all()
        Loads all required file inputs.
        Initiates values, called also by constructor.
    update(name, value, item=...)
        Update or add a data record with new values.
    delete(name, item)
        Delete given data item (data lists only).
    get_item(name, item=..)
        Return required data item.
    get_scan(name)
        Return scan data for use with :func:`stressfit.commands.set_scan`
    is_ready()
        Check if the input is ready for use (all file data loaded, 
        all associated data exist).

    """
    
    # all registered data sets
    reg_keys = ['shape','geometry','sampling','attenuation','data','imodel','emodel']
    # data lists
    list_keys = ['geometry','sampling','data','imodel','emodel']
    # items which require other loaded from files
    load_keys = ['sampling','attenuation','data']
    
    def __init__(self, config):
        # keys of listed data items
        self._keys = {}
        for key in UI_data.list_keys:
            self._keys[key] = []
        # flag for changed state
        self._reload_needed = {}
        # get initial listed data from resources
        self._data = config['udata']
        # get initial input data from resources
        self._input = config['uinput']        
        # check consistency
        self._validate_keys()
        self._update_keys()
        # load data from files
        self._force_reload()
        self.reload_all()


    def _validate_keys(self):
        for key in UI_data.reg_keys:
            if not key in self._data:
                raise Exception('Incomplete input data: {}.'.format(key))
        for key in self._data:
            if not key in UI_data.reg_keys:
                raise Exception('Undefined data group loaded: {}.'.format(key))
        for name in UI_data.list_keys: 
            kin = self._input[name].keys()
            for item in self._data[name]['list']:
                kdat = list(self._data[name]['list'][item].keys())
                if 'fdata' in kdat:
                    kdat.remove('fdata')
                for k in kdat:
                    if not k in kin:
                        raise Exception('Data has undefined item: {}[{}].{}.'.format(name, item, k))
                for k in kin:
                    if not k in kdat:
                        raise Exception('Data has missing item: {}[{}].{}.'.format(name, item, k))    

    def _update_keys(self):
        """Update references to list keys."""
        for key in UI_data.list_keys:
            self._keys[key].clear()
            self._keys[key].extend(list(self._data[key]['list'].keys())) 

    @property
    def inputs(self):
        """Return input data."""
        return self._input
    
    @property
    def data(self):
        """Return input data."""
        return self._data

    def _set_reload(self, name, value, item=None):
        """Set reload flag for given data item.
        
        Parameters
        ----------
        name : str
            Name of the data set.
        value : bool
            Flag value.
        item : str
            Name of the item for a listed data set.
            If None, set the flag to all items.
        """
        # skip items which do not require file input
        if not name in UI_data.load_keys:
            return
        # data is a list
        if name in UI_data.list_keys:
            if not name in self._reload_needed:
                self._reload_needed[name] = {}
            if item in self._data[name]['list']:
                self._reload_needed[name][item] = value
            else:
                self._reload_needed[name].clear()
                for key in self._data[name]['list']:
                    self._reload_needed[name][key] = value
        # data is a single item
        else:
           self._reload_needed[name] = value 

    def _force_reload(self):
        """Reset all reload flags to true.
        
        Forces reloading of all data on reload_all().
        """        
        self._reload_needed.clear()
        for k in UI_data.load_keys:
            self._set_reload(k, True)

    def _load_sampling(self, key):
        """Load sampling data."""
        try:
            sdata = self._data['sampling']['list'][key]
            sampling = comm.load_sampling(**sdata)
            sdata['fdata'] = sampling
        except Exception as e:
            sdata['fdata'] = None
            raise Exception(e)

    def _load_att(self):
        """Load attenuation data."""        
        if self._data['attenuation']['type'] == 'table':
            att = self._data['attenuation']['table']
            try:
                att['fdata'] = dataio.load_data(att['file'],
                                               path=att['path'],
                                               kind='tables',
                                               verbose=False)
            except Exception as e:
                att['fdata'] = None
                raise Exception(e)

    def _validate_scan(self, item):
        """Check input consistency for given scan.
        
        Parameters
        ----------
        item : str
            Key string for the data item.

        Returns
        -------
        bool, str[]
            Result of validation and error messages.
        """      
        msgs = []               
        dlist = self._data['data']['list']
        # named data set exist
        if not item in dlist:
            msgs.append('Input data item not defined: {}'.format(item))         
        # data item is present
        elif not 'fdata' in dlist[item]:
            msgs.append('No input data for "{}"'.format(item))
        # data is loaded
        elif not isinstance(dlist[item]['fdata'], ScanData):
            msgs.append('Input data not loaded for "{}"'.format(item))          
        # check association with geometry and sampling
        if item in dlist:
            ori = dlist[item]['geometry']            
            if not ori in self._data['geometry']['list']:
                msgs.append('No such geometry data: {}'.format(ori))
            else:
                gobj = self._data['geometry']['list'][ori]
                qry = [k in gobj for k in Geometry.input_keys]
                if not all(qry):
                    msgs.append('Geometry is not defined: {}'.format(ori))
            sam = dlist[item]['sampling']
            if not sam in self._data['sampling']['list']:
                msgs.append('No such sampling data: {}'.format(sam))
            else:
                sobj = self._data['sampling']['list'][sam]
                if sobj is None:
                    msgs.append('Sampling is not loaded: {}'.format(sam))
                    
        out = len(msgs)==0
        return out, msgs
    
    def _load_data(self, key):
        """Load input (measured) data.""" 
        vals = self._data['data']['list'][key]
        if vals['intensity']:
            ifile = vals['intensity'].strip()
        else:
            ifile = None
        sfile = vals['strain'].strip()  
        if not sfile or sfile=='.':
            vals['fdata'] = None
            msg = 'No strain file specified for [{}]'.format(key)
            raise Exception(msg) 
        try:
            scandata = ScanData()
            scandata.load(sfile, fint=ifile)
            vals['fdata'] = scandata
        except Exception as e:
            vals['fdata'] = None
            msg = 'Cannot load data [{}].'.format(key)
            msg += str(e)
            raise Exception(msg)          


    @property
    def items(self):
        """Return persistent lists of item keys for all data sets.
        
        Can be used as options for dropdown widgets. The list content is 
        up to date with actual data content. 
        
        """
        return self._keys
    
    @property
    def reload_needed(self):
        """Return true if input has changed and requires data reload."""
        out = False
        for name in self._reload_needed:
            if not out:
                if name in UI_data.list_keys:
                    qry = any(self._reload_needed[name].values())
                else:
                    qry = self._reload_needed[name]
                out = out or qry
        return out

    def input_from_dict(self,value):
        """Import input parameters from a dictionary."""
        if 'uinput' in value:
            lst = value['uinput']
            for key in self._input:
                if key in lst:
                    self._input[key] = lst[key]
        else:
            print('No uinput data found.')        

    def data_from_dict(self, value):
        """Import input data from a dictionary."""
        if 'udata' in value:
            lst = value['udata']
            for key in UI_data.reg_keys:
                if key in lst:
                    self._data[key] = lst[key]
            # check consistency
            self._validate_keys()
            self._update_keys()
        else:
            print('No uconfig data found.')

    def data_to_dict(self):
        """Export input data as dict."""       
        out = copy.deepcopy(self._data)
        # remove file data content
        for what in UI_data.list_keys:
            lst = out[what]['list']
            for item in lst:
                if 'fdata' in lst[item]:
                    lst[item]['fdata'] = None 
        out['attenuation']['table']['fdata'] = None
        return {'udata':out}

    def load(self, what, item=None):
        """Reload file input for given data set.
        
        Parameters
        ----------
        what : str
            Data set name (e.g. "sampling")
        item : str
            Name of the data list item. If None, reload all items of given list.
        
        """
        if not what in UI_data.load_keys:
            return
        try:
            if what == 'sampling':
                if item:
                    if item in self._data[what]['list']:
                        self._load_sampling(item)
                else:
                    for key in self._data[what]['list']:
                        self._load_sampling(key)                  
            elif what == 'attenuation':
                self._load_att()          
            elif what == 'data':
                if item:
                    if item in self._data[what]['list']:
                        self._load_data(item)
                else:
                    for key in self._data[what]['list']:
                        self._load_data(key)
            self._set_reload(what, False, item=item)
        except Exception as e:
            print(e)
    
    def clean(self):
        """Clean idle references etc.."""
        for what in UI_data.load_keys:
            if what in UI_data.list_keys:
                for key in list(self._reload_needed[what]):
                    if not key in self._data[what]['list']:
                        del self._reload_needed[what][key]
    
    def reload(self, name, item=None, force=False):
        """Reload required file data item."""
        if not name in UI_data.load_keys:
            return
        if name in UI_data.list_keys:
            if item in self._data[name]['list']:
                # add reload flag for given data item if missing
                if not item in self._reload_needed[name]:
                    self._reload_needed[name][item] = True
                if force or self._reload_needed[name][item]:
                    self.load(name, item=item)        
        else:
            if force or self._reload_needed[name]:
                self.load(name)
            
    def reload_all(self, force=False):
        """Load all required file data item."""
        for what in UI_data.load_keys:            
            if what in UI_data.list_keys:                
                for key in self._data[what]['list']:
                    self.reload(what, item=key, force=force)            
            else:
                self.reload(what, force=force)
        self.clean()

    def get_attenuation(self):
        """Return attenuation either as a table or scalar value.
        
        Compatible with :func:`stressfit.commands.set_attenuation`.

        """
        a = self._data['attenuation']
        if a['type'] == 'table':
            if not 'fdata' in a['table'] or a['table']['fdata'] is None:
                self._load_att()
            out = a['table']['fdata']
        else:
            out = a['value']
        return out
                
    def get_scan(self, name):
        """Return scan data as dict.
        
        Parameters
        ----------
        name : str
            Key string for the data item.

        Returns
        -------
        dict
            Scan data compatible with :func:`stressfit.commands.set_scan`.

        """
        valid, msgs = self._validate_scan(name)
        if valid:
            dset = self._data['data']['list'][name]
            sobj:ScanData = dset['fdata']
            ori = dset['geometry']
            sam = dset['sampling'] 
            geom = self._data['geometry']['list'][ori]
            sampling = self._data['sampling']['list'][sam]['fdata']
            sobj.assign_geometry(geom)
            sobj.assign_sampling(sampling)
            return sobj.scan
        else:
            for msg in msgs:
                print(msg)
            return None 
        
    def get_item(self, name, item=None):
        """Return requested input data.

        Parameters
        ----------
        name : str
            Name of the data set, one of UI_data.reg_keys.
        item : str, optional
            Key for the list item (only used with data lists)

        Returns
        -------
        dict
            Content of the requested data set. For data lists, return 
            the whole list if item=None.

        """
        if name in UI_data.list_keys:
            if item:
                return self._data[name]['list'][item]
            else:
                return self._data[name]['list']
        else:
            rec = self._data[name]
        
        return rec    

    def listkeys(self, key):
        """Return keys for given listed data.
        
        Use as a permanent link to a component which requires
        updated list of keys, e.g. options for a Dropdown list. 

        Parameters
        ----------
        key : str
            Key of corresponding data list, one of UI_data.list_keys.
        
        Returns
        -------
        list
            Keys for given data list.
        """
        if key in self._keys:
            out = self._keys[key]
        else:
            out = None
        return out
    

    def update(self, name, value, item=None):
        """Update input data or add new item to a list.
        
        For data lists, the content is updated with the value if
        the item already exists. Otherwise, a new item is added.

        Parameters
        ----------
        name : str
            Name of the data set, one of UI_data.reg_keys.
        value : dict
            Content of the input.
        item : str, optional
            Key for the list item (only used with data lists)

        """
        if name in UI_data.list_keys:
            lst = self._data[name]['list']
            if item is None:
                rec = lst
            elif item in lst:
                rec = lst[item]
            else:
                rec = {}
                lst[item] = rec
            rec.update(value)
            self._keys[name].clear()
            self._keys[name].extend(list(lst.keys()))
        elif name in UI_data.reg_keys:
            self._data[name].update(value)
        self._set_reload(name, True, item=item)
    
    def delete(self, name, item):
        """Delete specified item of the named data list."""
        lst = self._data[name]['list']
        if name in UI_data.list_keys and item in lst:
            del self._data[name]['list'][item]
            self._keys[name].clear()
            self._keys[name].extend(list(lst.keys()))
            if name in UI_data.load_keys and item in self._reload_needed[name]:
                del self._reload_needed[name][item]
    
    def validate(self):
        """Check input consistency for all input data."""
        msgs = []
        dlist = self._data['data']['list']
        valid = True
        for key in dlist:
            v, s = self._validate_scan(key)
            valid = valid and v
            msgs.extend(s)
        if not valid:
            for msg in msgs:
                print(msg)            
        return valid
  
    def is_ready(self):
        """Check input consistency and make sure file input is loaded."""
        out = True
        if not self.validate():
            out = False
            print('Fix the input.')
        elif self.reload_needed:
            out = False
            print('Input has changed, call reload_all() and check is_ready() again.')          
        return out
      

class UI_config():
    """Container for user configuration of stressfit.
    
    Contains all user input, user data and program configuration necessary to 
    run stressfit commands. 
    """
    
    def __init__(self):
        # get initial configuretion from resources
        self._config = dataio.load_config('config_ui')
        self._uconfig = self._config['uconfig']
        self._udata = UI_data(self._config)
        self._keys =  list(self._uconfig.keys())
    
    @property
    def keys(self):
        """Return allowed value keys."""
        return self._keys
    
    @property
    def data(self):
        """Return UI_data object."""
        return self._udata

    @property
    def config(self):
        """Return configuration parameters as dict."""
        return {'uconfig':self._uconfig}

    def export_dict(self):
        """Export all user input as a dictionary."""
        out = {}
        out.update(self.data.data_to_dict())
        out.update(self.config)
        out.update({'uinput':self.data.inputs})
        return out

    def export_json(self):
        """Export all user input in json format."""
        def clean(obj):
            """Compact lists."""
            s = ''
            L = len(obj.groups())
            if L>0:       
                g = obj.groups()[0]
                s = re.sub(r'\s{2,}',r' ',g)
            return s
        txt = self.export_dict()
        txt1 = json.dumps(txt,indent=4)
        out = re.sub(r'(\[[^]]*\])',clean,txt1)
        return out
        
    def import_dict(self, inp, reload=False):
        """Import all user input from a dictionary.
        
        Parameters
        ----------
        inp : dict
            User input as dictionary
        reload : boolean
            Reload required file data.
        """
        self.config_from_dict(inp) 
        self.data.input_from_dict(inp)
        self.data.data_from_dict(inp) 
        if reload:
            self.reload_all(force=True)

    def import_json(self, lines, reload=False):
        """Import all user input from a text in json format.
        
        Parameters
        ----------
        lines : list
            List of text lines.
        reload : boolean
            Reload required file data.
        """
        inp = json.loads('\n'.join(lines)) 
        self.import_dict(inp, reload=reload)

    def config_from_dict(self, value):
        """Import configuration parameters from a dictionary."""
        if 'uconfig' in value:
            lst = value['uconfig']
            for key in self._uconfig:
                if key in lst:
                    self._uconfig[key] = lst[key]
        else:
            print('No uconfig data found.')
    
    def get(self, key):
        """Return requested config data.

        Parameters
        ----------
        key : str
            Value name.

        Returns
        -------
        dict
            Content of the requested input.

        """
        rec = None
        if key in self._uconfig:
            rec = self._uconfig[key]        
        return rec
    
    def update(self, key, value):
        """Update input data.

        Parameters
        ----------
        key : str
            Value name.
        value : dict
            Content of the input.

        """
        if key in self._uconfig:
            self._uconfig[key].update(value)
    
    def validate(self, data:bool=True):
        """Check input consistency with UI data.
        
        Parameters
        ----------
        data: bool
            If true, validate also UI_data.
        """
        out = True
        if data:
            out = self._udata.validate()
            if not out:
                print('Fix config data input.')
                return out
        
        msgs = []
        # check that sampling and geometry values are define in the input data
        for key in self._uconfig:
            for opt in ['sampling','geometry']:
                if opt in self._uconfig[key]:
                    item = self._uconfig[key][opt]
                    dlist = self._udata.get_item(opt)
                    if not item in dlist:
                        out = False
                        msg = '{}: undefined {} [{}]'.format(key, opt, item)
                        msgs.append(msg)
        if not out:
            for msg in msgs:
                print(msg)         
        return out
    
    def is_ready(self):
        """Check that all inpt is valid and ready for use."""
        out1 = self._udata.is_ready()
        out2 = self.validate(data=False)
        return out1 and out2

    def reload_all(self, force=False):
        """Reload all file input."""
        self._udata.reload_all(force=force) 
        
    
#%%
def test():
    """Basic module test."""
    dataio.workspace().info()
    uiconf = UI_config()
    hoop = {
            "strain": "eps_SS_hoop.dat",
            "intensity": "int_SS_hoop.dat",
            "geometry": "radial",
            "sampling": "1x1x5mm",
    		"fdata": None
            }
    uiconf.data.update('data', hoop, item='hoop')
    res = uiconf.is_ready()
    print(res)
    
    uiconf.reload_all()
    ori = uiconf.data.get_item('geometry',item='radial').copy()
    ori['angles'] = [135,0,90]
    uiconf.data.update('geometry',ori,item='hoop')
    res = uiconf.is_ready()
    print(res)

