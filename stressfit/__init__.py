import numpy as np
# For each measured reflection and sample orientation, provide
# - strain data file: sample position (encoder value), strain, error
# - intensity data file: sample position (encoder value), intensity, error
# - scan direction and origin (encoder = 0) in sample coordinates
# - sample rotation centre (sample coordinates)
# - sample orientation (Euler angles)


def load_input(epsfile, intfile,
              path="",
              scandir=[0, 0, 1], 
              scanorig=[0, 0, 0], 
              rotctr=[0, 0, 0],
              angles=[180+45, 0, 0],
              sampling=None):
    """Load experimental data and metadata.
    
    Parameters
    ----------
    epsfile : str
        File name for strain data: position (encoder value), strain, error
    intfile : str, optional
        File name for intensity data: position (encoder value), strain, erro.
    scandir : list or array, optional
        Scan direction in sample coordinates
    scanorig : list or array, optional
        Sscan origin (encoder = 0) in sample coordinates
    rotctr : list or array, optional
        Sample rotation centre (sample coordinates)
    angles : list or array, optional
        Sample orientation (Euler angles YXY)
    sampling: dict
        Sampling events loaded by the function load_gauge.

    Returns
    -------
    dict
        Input data: keys are derived from argument names.

    """
    out = {}
    out['eps'] = np.loadtxt(path+epsfile)
    if intfile:
        out['int'] = np.loadtxt(path+intfile)
    else:
        out['int'] = None
    out['scandir'] = np.array(scandir)
    out['scanorig'] = np.array(scanorig)
    out['rotctr'] = np.array(rotctr)
    out['angles'] = np.array(angles)
    out['sampling'] = sampling
    return out


def load_sampling(filename, columns=[1, 4, 7, 10, 11], maxn=0, verbose=False):
    """Load Monte Carlo events representing the sampling distribution.
    
    The event table contains neutron coordinates, weights and dhkl values.  
    You need to specify column numbers for position, ki and kf vectors, 
    weights and dhkl.
    
    Imported MC events are defined in the laboratory frame, with the origin 
    at the centre of the instrumental gauge volume. Sample and orientation 
    thus defines zero scan position and scan direction in the sample.
    

    Parameters
    ----------
    filename : str
        Input file name.
    columns : list, optional
        Column indices of r[0], ki[0], kf[0], weight and dhkl (starts from 0)
    maxn : int, optional
        Maximum number or records. If zero, take all records in the file.
    verbose: boolean
        If true, print calculated gauge parameters.

    Returns
    -------
    dict
        - data: event data
        - columns: given input parameter
        - nrec: number of records
        - ctr: sampling centre
        - dmean: mean dhkl
    """
    out = {}
    data = np.loadtxt(filename)
    if maxn:
        nrec = min(maxn, data.shape[0])
    else:
        nrec = data.shape[0]
    # Calculate centre of mass of the distribution 
    P = data[:nrec,columns[3]]/np.sum(data[:nrec,columns[3]])
    ctr = np.zeros(3)
    ki = np.zeros(3)
    kf = np.zeros(3)
    for i in range(3):
        ctr[i] = data[:nrec, columns[0] + i].dot(P)
        ki[i] = data[:nrec, columns[1] + i].dot(P)
        kf[i] = data[:nrec, columns[2] + i].dot(P)
    
    wav = 2*np.pi/np.sqrt(ki.dot(ki))
    dmean = data[:nrec, columns[4]].dot(P)
    tth = np.arcsin(wav/2/dmean)*360/np.pi
    if verbose:
        print('Loaded event list with {:d} records'.format(nrec))    
        print('Gauge centre: [{:g}, {:g}, {:g}] '.format(*ctr))
        print('Mean wavelength: {:g}'.format(wav))
        print('2 theta: {:g}'.format(tth))
        print('d0: {:g}\n'.format(dmean))    
    out['data'] = data[:nrec,:]
    out['columns'] = columns
    out['nrec'] = nrec
    out['ctr'] = ctr
    out['dmean'] = dmean
    out['wav'] = wav
    out['tth'] = tth
    out['ki'] = ki
    out['kf'] = kf
    return out
