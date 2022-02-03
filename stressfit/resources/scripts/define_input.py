"""Script for setting stressfit input parameters.

Created on Tue Oct 26 16:31:49 2021
@author: Jan Saroun, saroun@ujf.cas.cz
"""
import stressfit.dataio as dataio
import stressfit.commands as comm

# define workspace
#--------------------------------------------
wks = dataio.workspace()
workspace = wks.get_paths()
# modify default workspace paths as needed
# 'work' directory must be absolute
# other paths can be relative to 'work' directory

# workspace['work'] = 'my new workspace direcotry'
# workspace['data'] = './input'
# workspace['tables'] = './tables'
# workspace['output'] = './out'

# apply changes
wks.set_paths(**workspace)

# sample shape
#--------------------------------------------
shape = {"shape": "Tube",
         "param": {"Rin": 5.0, 
                   "Rout": 9.0, 
                   "height": 50.0,
                   "ctr": [0.0, 0.0],
                   "sref": 0
                   }
         }

# list of geometries
#--------------------------------------------
geometry = {}
glist = {}
glist['radial'] = {"angles": [135, 0, 0],
                   "scandir": [0,0,-1],
                   "rotctr": [0,0,0],
                   "scanorig": [0,0,0]
                   }
geometry['list'] = glist

# list of sampling data
#--------------------------------------------
sampling = {}
slist = {}
slist['1x1x5mm'] = {"path": workspace['data'],
                    "file": "events_S_1mm.dat",
                    "nev": 3000
                    }
sampling['list'] = slist

# beam attenuation data
#--------------------------------------------
attenuation = {"table": {"path": workspace['tables'],"file": "Fe_mu.dat"},
               "value": 1.1,
               "type": "table"
               }

# collect
#--------------------------------------------
setup = {}
setup['workspace'] =workspace 
setup['shape'] = shape
setup['geometry'] = geometry
setup['sampling'] = sampling
setup['attenuation'] = attenuation
comm.set_user_input(setup)

