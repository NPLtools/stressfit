# -*- coding: utf-8 -*-
# Created on Wed Aug 25 16:05:34 2021
# @author: Jan Saroun, saroun@ujf.cas.cz
"""A module describing experiment geometry.

classes
-------
Geometry
    Encapsulates experiment geometry such as scan direction,
    sample rotation etc.

"""
import numpy as np
_deg = np.pi/180

class Geometry: 
    """Class describing experiment geometry.

    Describes properties such as sample orientation, 
    scan origin and direction etc ...
    
    The class assumes the basis frame for a usual neutron 
    diffractometer arrangement with ki || z and y-axis vertical.
    
    
    Attributes
    ----------
        scandir : array(3), optional
            Scan direction in sample coordinates. Internally, it defines
            moving of gauge across stationary sample.
            
            **NOTE** This is inverted with respect to the equivalent argument 
            of :func:`define`!
        scanorig : array(3), optional
            Scan origin (encoder = 0) in sample coordinates, in [mm]
        rotctr : array(3), optional
            Sample rotation centre (sample coordinates), in [mm]
        angles : array(3), optional
            Sample orientation (Euler angles YXY), in [rad]
        rot: array(3,3)
            Rotation matrix (transform vectors to sample frame)
    """
    
    input_keys = ['scandir', 'scanorig', 'angles', 'rotctr']
    
    def __init__(self, **kwargs):
        """Create Geometry object with default settings.
        
        See the :func:`~set` function for custom setting. 
        """
        self.define(**kwargs)
    
    def define(self, scandir=[0, 0, -1], scanorig=[0, 0, 0], angles=[0, 0, 0],
            rotctr=[0, 0, 0], **kwargs):
        """Define geometry properties.
        
        Parameters
        ----------
        scandir : array(3), optional
            Scan direction in sample coordinates (moving of events 
            across stationary sample)
        scanorig : array(3), optional
            Scan origin (encoder = 0) in sample coordinates, in [mm]
        angles : array(3), optional
            Sample orientation (Euler angles YXY), in [deg]
        rotctr : array(3), optional
            Sample rotation centre (sample coordinates), in [mm]
        """
        # NOTE: we must invert scandir, because it should define move of events
        # in stationary sample. 
        self.scandir = np.array(scandir)
        self.scanorig = np.array(scanorig)
        self.rotctr = np.array(rotctr)
        self.angles = np.array(angles)*_deg
        if 'name' in kwargs:
            self.name = kwargs['name']
        self.rot = Geometry.EulerYXY(self.angles)
        
    def toSample(self, v):
        """Transform a vector to the (rotated) sample frame."""
        return self.rot.dot(v-self.rotctr) + self.rotctr

    def fromSample(self, v):
        """Transform a vector from the (rotated) sample frame."""
        return self.rot.T.dot(v-self.rotctr) + self.rotctr  

    def EulerYXY(angles):
        """Return rotation matrix for given YXY Euler angles.
        
        Describes passive rotation: r' = R.r 
        where r is in the basis coordinates and r' is the same vector
        expressed in the rotated basis. 
        
        Parameters
        ----------
        angles: array(3)
            YXY Euler angles [rad]
        
        """
        R1 = Geometry.rmat(angles[0],axis=1)
        R2 = Geometry.rmat(angles[1],axis=0)
        R3 = Geometry.rmat(angles[2],axis=1)
        R = R3.dot(R2.dot(R1)) 
        return R


    def rmat(angle, axis=1):
        """Return elementary rotation matrix for given basis direction.
        
        Describes passive rotation: r' = R.r 
        where r is in the basis coordinates and r' is the same vector
        expressed in the rotated basis. 
        
        Parameters
        ----------
        angle: float
            rotation angle [rad]
        axis: int
            basis axis index (x=0, y=1, z=2)
        
        """
        s, c = np.sin(angle), np.cos(angle)
        if (axis == 1):
            R = np.array([[c, 0., -s], [0., 1., 0.], [s, 0., c]])
        elif (axis == 0):
            R = np.array([[1., 0., 0.], [0., c, s], [0., -s, c]])
        elif (axis == 2):
            R = np.array([[c, s, 0.], [-s, c, 0.], [0., 0., 1.]])
        else:
            R = np.eye(3)
        return R
    

    def rotate(v, axis, angle):
        """Rotate a vector around given axis.
    
        Parameters
        ----------
            v : array(3)
                3D vector to be rotated
            axis: int
                rotation axis index 0 .. 2
            angle: float
                rotation angle in rad
        Returns
        -------
            rotated vector
        """
        R = Geometry.rmat(angle, axis=axis)
        vv = np.array(v).reshape((3, 1))
        r = R.T.dot(vv).reshape(3,)
        return r
        
