# -*- coding: utf-8 -*-
# Created on Wed Aug 18 21:09:03 2021
# @author: Jan Saroun, saroun@ujf.cas.cz
"""STRESSFIT module which deals with lattice response on applied stress."""

import numpy as np
import stressfit.dataio as dataio
import json

class Compliance():
    """Class representing compliance of material.
    
    Compliance is internally represented by the S_33kl subset of the 
    compliance tensor in the Q-coordinates (Q||z), plus corresponding rotation 
    matrix R transforming vectors from Q-frame to sample frame.
    
    Compliance is hkl sensitive due to averaging over all grain orientations 
    with Q || [hkl].
    
    Provides stress factors F_i which are used to calculate measured lattice 
    strain from stress tensor in sample reference frame: 
        eps = F_i*sigma_i
    where i=0..5 is the tensor indexing in Voigt notation. 
    
    In STRESSFIT, the sample frame is aligned with the sample shape 
    (see documentation for individual Shape objects). 
    This class neglects possible misalignment between the sample frame and 
    principal stress directions. This means that the input stress tensor 
    is always assumed to be expressed in the sample reference frame.
    
    Attributes
    ----------
    F: array(6,:)
        Stress factors in Voigt notation, for calculation of measured strain 
        from stress tensor in sample reference frame. Units: 1e-6/MPa
    R: list of array(3,3)
        Rotation matrix, transforming vectors from Q-reference (Q||z) 
        to sample reference frame.
    hkl: list of array(3)
        hkl indices for the r.l. direction corresponding to the stress factors.
    Q: list of array(3)
        Normalized Q-vector in sample coordinates.
    texture: string
        Information about texture.
    method: string
        Information about averaging method.

    """
    
# static methods  
  
    def normalize(v):
        """Return normalized vector v."""
        a = np.reshape(v,(-1,3))
        v0 = np.sqrt(np.sum(a*a, axis=1))
        out = np.divide(a,v0.reshape((-1,1)))
        return out
    
    def S33toV6(m):
        """Convert S_33kl compliance matrix to 6-vector in Voigt notation.
        
        Includes the multiplication rule: V(i) = 2*M(k,l) for k<>l.
        i=3 for k,l=1,2 etc.
        """
        V = np.zeros(6)
        iV = [[0,0],[1,1],[2,2],[1, 2],[0,2],[0,1]]
        f = [1.0,1.0,1.0,2.0,2.0,2.0]
        for i in range(6):
            V[i] = f[i]*m[iV[i][0],iV[i][1]]
        return V
    
    def V6toS33(V):
        """Convert 6-vector in Voigt notation to S_33kl compliance matrix.
        
        Includes the multiplication rule: V(i) = 2*M(k,l) for k<>l.
        i=3 for k,l=1,2 etc.
        """        
        m = np.zeros((3,3))
        iV = [[0,0],[1,1],[2,2],[1, 2],[0,2],[0,1]]
        f = [1.0,1.0,1.0,2.0,2.0,2.0]
        for i in range(6):
            k=iV[i][0]
            l=iV[i][1]
            m[k,l] = V[i]/f[i]
            m[l,k] = m[k,l]
        return m
    
    def s_from_v6(v, rot=None):
        """Stress tensor from Voigt 6-vector."""
        iV = [[0,0],[1,1],[2,2],[1, 2],[0,2],[0,1]]
        m = np.zeros((3,3))
        for i in range(6):
            k=iV[i][0]
            l=iV[i][1]
            m[k,l] = v[i]
            m[l,k] = m[k,l]
        if rot is not None:
            out = rot.dot(m.dot(rot.T))
        else:
            out = m
        return out
    
    def v6_from_s(s, rot=None):
        """Voigt 6-vector from stress tensor."""
        V = np.zeros(6)
        iV = [[0,0],[1,1],[2,2],[1, 2],[0,2],[0,1]]
        for i in range(6):
            V[i] = s[iV[i][0],iV[i][1]]
        return V    
    
    def angles_to_q(angles):
        """Calculate Q-vectors from orientation angles phi, psi.
        
        Parameters
        ----------
        angles: list or array(2)
            Input angles [phi, psi] in rad. 
        
        Returns
        -------
        list of Q-vectors
            Normalized directions of the scattering vector. 
        """
        a = np.array(angles).reshape((-1,2))
        s = np.sin(a)
        c = np.cos(a)
        qx = c[:,0]*s[:,1]
        qy = s[:,0]*s[:,1]
        qz = c[:,1]
        out = np.array([qx,qy,qz]).T
        return out
    
    
    def q_to_angles(Q):
        """Calculate orientation angles phi, psi from Q vector.
        
        Parameters
        ----------
        Q: list or array(3)
            Input Q-vectors given as array[:,3] or list of array(3)
        
        Returns
        -------
        array of [phi, psi]
            ZYZ Euler angles in [rad] (set the 3rd angle to 0).
            Returns a pair of angles for each Q-vector.
        """
        qn = Compliance.normalize(Q)
        phi = np.arctan2(qn[:,1], qn[:,0])
        psi = np.arccos(qn[:,2])
        out = np.array([phi,psi]).T
        return out

# constructor and class methods

    def __init__(self, hkl, S33kl, angles):
        """Create instance of Compliance object.
        
        To be used in general case. See also other constructors:
            createIsotropic()
                For isotropic materials, using diffraction elastic
                constants as Yound modulus and Poisson constant.
            fromFile()
                Loads the input parameters from a file.
            :func:`load_resource`
            :func:`set_hkl
        
        Parameters
        ----------
        hkl: array(3)
            hkl values for corresponding reflection.
        S33kl: array(6)
            S_33kl elements of the compliance tensor in Q-coordinates (Q|z).
        angles: array(2)
            ZYZ Euler angles in [rad]. They define ZYZ rotation matrix which 
            transforms vectors from Q-coordinates to sample frame. Hence
            Q = R.[0,0,1].T. The 3rd angle is irrelevant, hence it is set to 
            zero so that y remains in the sample xy plane. Therefore, only 
            2 elements [phi, psi] are given on input.

        """
        self.hkl = hkl
        self.S33kl = S33kl
        self.set_q_angles(angles)
        self.texture = 'none'
        self.method = 'none'
        self.dbase = None


    @classmethod
    def create_isotropic(self, E_hkl, nu_hkl, hkl, Q=[0,0,1]):
        """Create Compliance instance for isotropic material.
        
        Parameters
        ----------
        E_hkl : float
            Young modulus, GPa
        nu_hkl : float
            Poisson constant
        hkl: array(3)
            hkl values for corresponding reflection.
        Q: array(3)
            Scattering vector direction in sample reference frame.
            Normalization is done by this function. Q can be later overriden
            by calling set_q() or set_q_angles(). 

        Returns
        -------
        Instance of Compliance for isotropic material.
        
        Note
        ----
        Compliance does not depend on Q for isotropic materials, but Q is 
        used to calculate measured strain from stress tensor in sample 
        coordinates.

        """
        # 3rd row of the compliance tensor in Q-reference, Voigt notation
        # Units: 1e-6/MPa
        V = 1000*np.array([-nu_hkl, -nu_hkl, 1, 0, 0, 0])/E_hkl
        # convert to tensor S_33kl
        S33_iso = Compliance.V6toS33(V)
        # calculate angles from Q-vector
        angles = Compliance.q_to_angles(Q)
        obj = self(hkl, S33_iso, angles)
        obj.texture = 'uniform'
        obj.method = 'none'
        return obj
    

    @classmethod
    def from_resource(self, filename, phase='', **kwargs):
        """Create compliance tensor object from resource file.
        
        Load the resource file and set select the compliance data
        for the first record found for the given phase. If phase is empty,
            select the first available phase record.

        Parameters
        ----------
        filename : str
            Resource file name (json format). See :func:`load_resource`
        phase : str, optional
            Phase name to be selected from the resource file. 
        **kwargs : 
            Optional arguments passed to :func:`load_resource`


        """
        obj = self(200, np.eye(3,3), [0,0])
        obj.load_resource(filename, **kwargs)
        item = obj._dbase_pick_first(phase)
        obj.set_hkl(item['phase'], item['hkl'], Q_angles=item['Q_angles'])
        return obj    

    @classmethod
    def from_file(self, filename, **kwargs):
        """Create compliance tensor object from file.
        
        Parameters
        ----------
        filename : string
            file name (full path)
        kwargs:
            arguments passed to load_from_file()

        """
        obj = self(200, np.eye(3,3), [0,0])
        obj.load_from_file(filename, **kwargs)
        return obj
    
    
# private methods
    def _dbase_info(self):
        """Return ifo string about current dbase content."""
        if self.dbase is None:
            return 'Database not defined.'
        rs = self.dbase
        hdrs = ['texture', 'method']
        out = 'filename: {}\n'.format(rs['filename'])
        for key in rs['phases']:
            val = rs['phases'][key]
            out += 'phase: {}\n'.format(key)
            for h in hdrs:
                if h in val:
                    out += '\t{}: {}\n'.format(h, val[h])
            fmt = '\t\thkl=[{:g},{:g},{:g}], Q=[{:g},{:g},{:g}]'
            fmt += ', Q_angles=[{:g},{:g}]\n'
            for item in val['data']:
                out += fmt.format(*item['hkl'],*item['Q'],*item['Q_angles'])
        return out

    def _dbase_pick(self, phase, hkl, **kwargs):
        out = None
        if self.dbase is None:
            raise Exception('Database not defined.')
        rs = self.dbase
        if not phase in rs['phases']:
            raise Exception('{} not found in {}'.format(phase, self.dbase['filename']))
        val = rs['phases'][phase]
        
        if 'Q_angles' in kwargs:
            Q_angles = kwargs['Q_angles']
        else:
            if 'Q' in kwargs:
                Q =  kwargs['Q']
            else:
                Q = [0,0,1]
            Q_angles = list(Compliance.q_to_angles(Q)*180/np.pi)

        for item in val['data']:
            if out is None:
                dif = np.abs(np.array(hkl) - np.array(item['hkl']))
                qry1 = dif<1e-2
                dif = np.abs(np.array(Q_angles) - np.array(item['Q_angles']))
                qry2 = dif<1e-2
                if qry1.all() and qry2.all:
                    out = item.copy()
                    out['Q_angles'] = Q_angles
                    # copy other properties to the output
                    args = ['texture','method']
                    for a in args:
                        if a in val:
                            out[a] = val[a]
                        
        if out is None:
            msg = 'Could not find requested record in {}.'
            raise Exception(msg.format(self.dbase['filename']))
        
        
        return out
    
    def _dbase_pick_first(self, phase=''):
        out = None
        if self.dbase is None:
            raise Exception('Database not defined.')
        rs = self.dbase
        if not phase in rs['phases']:
            p = list(rs['phases'].keys())[0]
            if phase:
                msg = 'WARNING: phase {} not defined. Using {}.'
                print(msg.format(phase,p))
            phase = p
        p = rs['phases'][phase]
        item = p['data'][0]
        item['phase'] = phase
        args = ['texture','method']
        for a in args:
            if a in p:
                item[a] = p[a]
        return item


# other methods   

    def info(self):
        """Print information on the current object setting."""
        print(self._dbase_info())


    def set_hkl(self, phase, hkl, **kwargs):
        """Set material data for new phase, hkl and Q from resources.
        
        Given phase and hkl and Q must be found in the previously loaded
        resource file (see :func:`load_resource`).
        
        Either Q_angles or Q should be defined. If both are defined, 
        than Q_angles are used.

        Parameters
        ----------
        phase: str
            Key for required phase from the resources.
        hkl : array(3)
            Required hkl values. Note that the index order must match
            the one given in the resource file, equivalent planes are not 
            recognized.
        
        kwargs
        ------
        Q_angles : array(2), optional
            Q orientation angles in [deg]. The values must match the input in 
            the resource file to 0.01 deg precission. Default is [0,0].
            See also :func:`~.set_q_angles`.
        Q : array(3), optional
            Required Q vector direction. The values must match the input in 
            the resource file to 1e-3 precission. Default is [0,0,1].
            If Q_angles is also defined, then Q is not used.

        Example
        -------
            `load_resource('Fe_iso')`
            
            `set_hkl('ferrite', [2,1,1], Q_angles=[0, 90])`

        """
        try:
            item = self._dbase_pick(phase, hkl, **kwargs)
            self.hkl = item['hkl']
            self.S33kl = Compliance.V6toS33(item['S'])
            self.set_q_angles(item['Q_angles'], unit='deg')
            if 'texture' in item:
                self.texture = item['texture']
            if 'method' in item:
                self.method = item['method']
        except Exception as e:
            print(e)
        
        
    def set_q_angles(self, ori, unit='rad'):
        """
        Set new orientations of Q-vector in sample frame.

        Parameters
        ----------
        ori : list of array(2)
            Phi and psi form the set of ZYZ Euler angles [Phi, psi, 0]. 
            Corresponding ZYZ rotation matrix (R) transforms vectors 
            from Q-coordinates to sample frame. It defines orientation 
            of the scattering vector in sample frame: Q = R.[0,0,1].T. 
            The 3rd Euler angle is irrelevant and is set to zero so that 
            y remains in the sample xy plane. 
        unit: str
            Set angular unit to 'rad' or 'deg'

        """
        # ori can be passed either as a single pair [phi,psi] 
        # or a list of [psi, phi] 
        # make sure we take it as a list of array(2) 
        uni = 1
        if unit=='deg':
            uni = np.pi/180             
        o = np.reshape(ori,(-1,2))
        nr = np.shape(o)[0]
        self.R = nr*[0]
        self.Q = nr*[0]
        self.F = np.zeros((6,nr))
        for i in range(nr):
            a = o[i,:]*uni
            s = np.sin(a)
            c = np.cos(a)
            self.R[i] = np.array([
                     [c[0]*c[1] , -s[0] , c[0]*s[1]],
                     [s[0]*c[1] ,  c[0] , s[0]*s[1]],
                     [    -s[1] ,  0    , c[1]]
                     ])
            # F = R.S.R^T
            F = self.R[i].dot(self.S33kl.dot(self.R[i].T))
            # store as Voigt 6-vector
            self.F[:,i] = Compliance.S33toV6(F)
            # Q = R.[0,0,1]
            self.Q[i] = self.R[i].dot(np.array([0,0,1]))
    
    
    def set_q(self, q):
        """Set Q orientations.
        
        Parameters
        ----------
        q : array
            Can be a single q=[qx, qy, qz] vector or a list of q.

        """
        ori = Compliance.q_to_angles(q)
        self.set_q_angles(ori)
    

    def load_resource(self, filename, path=None, verbose=True):
        """Load material compliance data from resources in JSON format.
        
        See resources/compliance for examples. Specify either full path,
        or a file name from package resources. Alternative path resources path
        can be provided as additional argument. 
        (See package resources/compliance for examples and file format.)

        Parameters
        ----------
        filename : str
            File name of the json database. Se
        path: str, optional
            Path to search for filename if relative
        verbose : boolean, optional
            Print content summary after loading.


        """
        def check_dbase(rs):
            todeg = 180/np.pi
            for p in rs['phases']:
                phase = rs['phases'][p]
                if not 'data' in phase:
                    msg = 'Wrong file format, `data` key not found for: {}'
                    raise Exception(msg.format(phase))
                for item in phase['data']:  
                    #item = phase['data'][key]
                    if (not 'Q_angles' in item) and (not 'Q' in item):
                        Q_angles = [0,0]
                        Q = [0,0,1]
                    elif not 'Q_angles' in item:
                        Q = item['Q']
                        Q_angles = Compliance.q_to_angles(Q)*todeg
                        
                    elif not 'Q' in item:
                        Q_angles = item['Q_angles']
                        Q = Compliance.angles_to_q(np.array(Q_angles)*todeg)
                    else:
                        Q_angles = item['Q_angles']
                        Q = item['Q']
                    item['Q_angles'] = list(Q_angles) 
                    item['Q'] = list(Q)
            
        if not (filename.endswith('.json') or filename.endswith('.JSON')):
            filename += '.json'
        s = dataio.load_text(filename, kind='tables', path=path)
        rs = json.loads('\n'.join(s))
        rs['filename'] = filename
        check_dbase(rs)
        self.dbase = rs
        if verbose:
            print(self._dbase_info())
            
        

    def load_from_file(self, filename, hkl=None, q=None, verbose=True):
        """Load compliance tensor from a text file.
        
        The file format assumes headers starting with #, 
        followed by 6x6 matrix. The matrix values should be in units of 1e-6/MPa.  
        
        Properties can be passed in the header in the format '# name: value'.
        The recogized properties are:
        
            hkl: h,k,l
                Miller indices for corresponding reflection
            q: (array of) [Qx, Qy, Qz]
                Q-vector in sample coordinates, can be an array of Qs
            Texture: string
                Information about texture (optional)
            Method: string
                Information about averaging method (optional)
        
        If `q` or `hkl` is passed as an argument,it overrides the file content.
        
        Parameters
        ----------
        filename : string
            file name (full path or a name from package resources)
        hkl:
            hkl indices (optional)
        q: (array of) [Qx, Qy, Qz]
                Q-vector in sample coordinates, can be an array of Q-vectors
        verbose: boolean
            print information about loading process

        """
        def str2v3(s):
            out = np.zeros(3)
            s = s.split(',')
            if len(s)==3:
                out = np.array(s).astype(np.float)
            else:
                raise Exception('Wrong format of vector string, expected x,y,z.')
            return out
                
        
        keys = ['hkl','z','texture','method']
        
        lines = dataio.load_text(filename, kind='tables')
        irow = 0
        S = np.zeros(6)
        my_hkl = hkl
        my_q = q
        for L in lines:
            LL = L.strip()
            # handle comments
            if LL.startswith('#'):
                s = LL[1:].split(':')
                if len(s)>1:
                    nam = s[0].strip().lower()
                    val = s[1].strip()
                    if nam in keys:
                        if nam=='hkl':
                            if not my_hkl:
                                my_hkl = str2v3(val)
                        elif nam=='z':
                            if not my_q:
                                my_q = str2v3(val)   
                        elif nam=='texture':
                            self.texture = val
                        elif nam=='method':
                            self.method = val
                        if verbose:
                            print('{}= {}'.format(nam, val));
            elif len(LL)>0 and irow<1:
                srow = LL.split();
                if len(srow)==6:
                    row = np.array(srow).astype(np.float)
                    S = row
                    irow += 1
        if irow != 1:
            raise Exception('Wrong format, expected one row of 6 floats.')    
        
        qry = (my_hkl is not None) and (my_q is not None)
        if qry:
            self.hkl = my_hkl
            self.S33kl = Compliance.V6toS33(S)
            self.set_q(my_q)
        else:
            raise Exception('Missing hkl or q vector in input.')
        if verbose:
            print('S matrix loaded from {}.'.format(filename))


    def strain(self, sigma):
        """Calculate strain tensor from given stress - fast version.
        
        Parameters
        ----------
        sigma: array(:,6)
            Stress tensors in MPa as 6-vectors (Voigt notation).
        
        Returns
        -------
        array
            Latice strains in units of 1e-6.
            The return array contains:
                - One row for each sigma 
                - One column for each previously defined orientation
                  (see `set_q` and `set_q_angles`).
        
        """
        ns = sigma.shape[1]
        out = sigma.dot(self.F[:ns,:])
        return out
        
    
    def strain_at_Q(self, sigma, Q=None):
        """Calculate strain tensor from given stress - general version.

        Parameters
        ----------
        sigma: array()
            Stress tensor in MPa, either as full 6-vector, or just 3-vector 
            with principal components.
        Q : array(:,3), optional
            New Q-directions.

       
        Returns
        -------
        array
            Latice strain values in units of 1e-6. 

        """
        if Q is not None:
            self.set_q(Q)
        sh = np.shape(sigma)
        nc = sh[len(sh)-1] 
        # make sure sig is 2D array with 6 columns
        sig = np.array(sigma).reshape((-1,nc))
        if nc==3:
            sig = np.concatenate((sig,np.zeros(sig.shape)),axis=1)
        eps = self.strain(sig) 
        sz = np.size(eps)
        if np.size(eps)==1:
            out = eps[0,0]
        elif sz==len(self.Q):
            out = eps[0,:]
        else:
            out = eps
        return out 
        
   
"""
class LatticeResponse():
    def __init__(self, sigmaL, sigmaT):
        self.sigmaL = sigmaL
        self.sigmaT = sigmaT
        intpL = interp1d(sigmaL[:,1], sigmaL[:,0], kind='linear')
        intpT = interp1d(sigmaT[:,1], sigmaT[:,0], kind='linear')
        sigmin = max(min(sigmaL[:,1]), min(sigmaT[:,1]))
        sigmax = min(max(sigmaL[:,1]), max(sigmaT[:,1]))
        sig = np.linspace(sigmin, sigmax, num=50)
        epsL = intpL(sig)
        epsT = intpT(sig)
        self.intpEL = interp1d(sig, epsL, kind='linear')
        self.intpET = interp1d(sig, epsT, kind='linear')
    
        
    def strain(self, sigma):
        s = np.reshape(sigma, (-1, 3))
        EL = self.intpEL(s)
        ET = self.intpET(s)
        exx = EL[:,0] + ET[:,1] + ET[:,2]
        eyy = EL[:,1] + ET[:,0] + ET[:,2] 
        ezz = EL[:,2] + ET[:,0] + ET[:,1]
        return np.array([exx,eyy,ezz]).T
"""

#%% Tests

def test():
    Q = [1,1,0]    
    C = Compliance.create_isotropic(200, 0.33, [1,1,1], Q=Q)
    #print('Q:\n',C.Q)
    #print('S33kl:\n',C.S33kl)
    
    # test 1
    s = np.array([0.5, 0.5, 0, 0, 0, 0.5])
    eps = C.strain_at_Q(s)
    ans = C.S33kl[2,2]
    cond = abs((eps-ans)/eps)<1e-6
    assert cond, 'failed test 1: \n{}\n{}'. format(eps, ans)  
    
    # test 2 
    C.set_q_angles([[45, 90],[135, 90], [0, 0]], unit='deg')
    s1 = [0, 0, 1]
    s2 = [0, 0, 2]
    s = np.array([s1,s2])
    eps = C.strain(s)
    ans = np.array([np.diag(C.S33kl),2*np.diag(C.S33kl)])
    cond = np.abs(np.subtract(eps, ans))<1e-8
    assert cond.all(), 'failed test 2: \n{}\n{}'. format(eps, ans)
    
     # test 3 
    C = Compliance.from_resource('compliance/Fe_duplex_1.4462_2018_20mm.json',
                                 phase='ferrite')
    # C.load_resource('compliance/Fe_duplex_1.4462_2018_20mm.json')
    C.set_hkl('austenite', [3,1,1], Q_angles=[0, 90])
    C.set_q_angles([0, 90], unit='deg')
    s1 = [0, 0, 1]
    s2 = [1, 0, 0]
    s = np.array([s1,s2])
    eps = C.strain(s)    
    ans = np.array([[C.S33kl[0,0]],[C.S33kl[2,2]]])
    cond = np.abs(np.subtract(eps, ans))<1e-8
    assert cond.all(), 'failed test 3: \n{}\n{}'. format(eps, ans)

