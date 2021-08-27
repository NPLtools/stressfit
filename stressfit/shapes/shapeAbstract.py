# -*- coding: utf-8 -*-
"""
Base abstract class for sample shape definition
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from numpy.linalg import norm
from . import tracing as tr
# from numpy.linalg import norm


class ShapeAbstract:
    shape_type = 'abstract'

    def __init__(self):
        self._R = np.eye(3, 3)
        self._isRot = False
        self._pos = np.zeros(3)

# Abstract methods to be overriden

    def depthLocal(self, rloc):
        """Return depths under the surfaces for given array of positions.

        Arguments:
            rloc -- array: positions in local coordinates
        Returns:
            [depth1, depth2, inside], where
            depth1 -- array: depth under the front surface (shape dependent)
            depth2 -- array: depth under the rear surface (shape dependent)
            inside -- array: 0 or 1 if outside or inside
        """
        raise NotImplementedError("Subclass must implement abstract method")

    def cross(self, r, k):
        """Calculate times to cross-sections with boundaries.

        Arguments:
            r -- position coordinates
            k -- direction vector
        Returns:
            List of [[t11, t12], [t21, t22], ...] arrays, one item per position
            [t11, t22] are entry and exit times, one row per layer encountered
            The times must be sorted (t11 < t12 < t21 < t22 ...)
        """
        raise NotImplementedError("Subclass must implement abstract method")

    def rayPaths(self, r, ki, kf):
        """Calculate input and output paths for given array of points.

        Arguments:
            r -- position coordinates
            ki -- input wave vector
            kf -- output wave vector
        Returns:
            [path1, path2]
        """
        raise NotImplementedError("Subclass must implement abstract method")

# Class methods which should work for all descendants

    def depth(self, r):
        """Return depths under surface in lab coordinates"""
        v = self.getLocalPos(r)
        return self.depthLocal(v)

    def scanDepth(self, x, xdir):
        """Return depth values for given scan positions and direction"""
        r = (x*np.array(xdir).reshape((3, 1))).T
        rloc = self.getLocalPos([0., 0., 0.])
        r1 = r + rloc
        [d, d2, ins] = self.depthLocal(r1)
        return d

    def pathInOut(self, r, ki, kf):
        """Calculate in and out paths through the shape

        Arguments:
            r -- position coordinates
            ki -- incident direction
            kf -- final direction
        Returns:
            list [ic, pin, pout]
            ic -- True if the position is inside the shape
            pin -- total input path through the shape interior
            pout -- total output path through the shape interior
        """
        # r1 = self.getLocalPos(r)
        # ki1 = self.getLocalDir(ki)/norm(ki)
        # kf1 = self.getLocalDir(kf)/norm(kf)
        cmi = self.cross(r, ki/norm(ki))
        cmf = self.cross(r, kf/norm(kf))
        # cmi and cmf have the same number of items equal to r.shape[0]
        rang = range(r.shape[0])
        res = []
        for i in rang:
            path1 = tr.crossPaths(cmi[i])
            is1 = path1[0]
            pin = path1[2]
            path2 = tr.crossPaths(cmf[i])
            is2 = path2[0]
            pout = path2[3]
            # TODO: Cases ic1 != ic2 should not happen.
            # Raise error?
            ic = (is1) and (is2)
            res.append([ic, pin, pout])
        return res

    def getRotation(self, omega, chi, phi):
        """ Generate YXY Euler matrix """
        s, c = np.sin(omega), np.cos(omega)
        R1 = np.array([[c, 0., -s], [0., 1., 0.], [s, 0., c]])
        s, c = np.sin(chi), np.cos(chi)
        R2 = np.array([[1., 0., 0.], [0., c, s], [0., -s, c]])
        s, c = np.sin(phi), np.cos(phi)
        R3 = np.array([[c, 0., -s], [0., 1., 0.], [s, 0., c]])
        R = R3.dot(R2.dot(R1))
        return R

    def reset(self):
        """Set position and orientation to zero"""
        self._R = np.eye(3, 3)
        self._isRot = False
        self._pos = np.zeros(3)

    def rotate(self, omega, chi, phi):
        """ Rotate shape using YXY Euler system """
        R = self.getRotation(omega, chi, phi)
        self._R = R.dot(self._R)
        self._pos = self._R.dot(self._pos)
        self._isRot = self._R.trace() < 3.0

    def moveTo(self, r):
        """Move shape by r """
        self._pos += np.array(r)
        
    def moveToAbs(self, r):
        """Move shape to r """
        self._pos = np.array(r)

    def getLocalPos(self, v):
        """ Transform position vector to local coordinates

        Returns:
            array(3, )
        """
        a = self.getLocalDir(v) - self._pos
        return a

    def getLocalDir(self, v):
        """ Transform direction vector to local coordinates

        Returns:
            array(3,)
        """
        if (self._isRot):    # rotate only if necessary
            r = self._R.dot(np.array(v))
        else:
            r = v
        return r
