# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""Class encapsulatig stressfit configuration and input data.

Implements class Config with following features:
    
    - manage all inputs for stressfit needed to execute stressfit functions
    - manage user workspace (input/output folders etc.]
    - maintain cross-dependences
    - validate input data 
    
"""
import stressfit.dataio as dataio
import json

class Config():
    
    # define allowed setup keys
    _keys = ['options',
             'shape',
             'geometry',
             'sampling',
             'attenuation',
             'scene',
             'resolution',
             'data']
    
    def __init__(self):
        self.wks = dataio.workspace()
        self._paths = {}  # workspace paths
        self._setup = {}  # program setup parameters (user input)
        self._data = {}   # input data (loaded from resources, calculated etc.)
        self.update_paths()

    def set_path(self, paths):
        self._paths.update(paths)
        
    def get_path(self, keys=None):
        if keys is None:
            return self._paths
        elif isinstance(keys, list):
            return {key:self._paths[key] for key in keys}
        else:
            return self._paths[keys]

    def set_setup(self, setup):
        self._setup.update(setup)
        
    def get_setup(self, keys=None):
        if keys is None:
            return self._setup
        elif isinstance(keys, list):
            return {key:self._setup[key] for key in keys}
        else:
            return self._setup[keys]

    def load(self, filename=''):
        """Load input data in JSON format."""

        if filename:
            fn = filename
        else:
            fn = self.wks.inputfile
        with open(fn, 'r') as f:
            lines=f.readlines()
        inp = json.loads('\n'.join(lines))
        if not 'stressfit' in inp:
            msg = 'Invalid stressfit configuration file: {}.'
            raise Exception(msg.format(fn))
        stp = inp['stressfit']
        if 'workspace' in stp:
            self._paths.update(stp['workspace'])
            self.validate_workspace()
        if 'setup' in stp:
            self._setup.update(stp['setup'])

            
    def save(self, filename='', workspace=False):
        """Save input data in JSON format."""
        if filename:
            fn = filename
        else:
            fn = self.wks.inputfile  
        out = {'setup': self._setup}
        if workspace:
            out['workspace'] = self._paths
        txt = json.dumps({'stressfit':out},indent=4)
        with open(fn, 'w') as f:
            f.write(txt)

            
    def update_paths(self):
        """Update workspace paths from dataio.workspace()."""
        self._paths.update(self.wks.get_paths())
    
    def update_workspace(self):
        """Update dataio.workspace according to config values."""
        wks_path = self.wks.get_paths('work')
        if wks_path != self._paths['work']:
            self.wks.change_workspace(self._paths['work'])
        else:
            self.wks.set_paths(**self._paths)
                
    def validate_workspace(self, verbose=False):
        """Check workspace configuration (paths must exist)."""
        try:
            self.update_workspace()
            self.wks.validate_paths()
            if verbose: 
                print('Workspace OK')
        except Exception as e:
            print(e)

    def reset_workspace(self):
        """Set workspace configuration to package default."""
        self.wks.reset_paths()
        self.update_paths()

    def save_workspace(self, default=False):
        """Save workspace information."""
        self.wks.save(default=default)

    def init(keys):
        
        if 'sampling' in keys:
            sampling = self._data[''] 
            comm.set_sampling(self._s['sampling'])
            if 'geometry' in kwargs:
                comm.set_geometry(kwargs['geometry'])