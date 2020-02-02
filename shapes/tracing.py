# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:12:32 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
# from numpy.linalg import norm
_eps = 1.0e-12
_inf = 1.0e+12


def crossLayer(w, r, k, ix):
    """Calculate cross-times with plane parallel layer

    Arguments:
        w -- layer thickness
        r -- array[:,3] of position vectors
        k -- direction vector
        ix --  dimension index [0 .. 2]
    Returns:
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
    """
    if abs(k[ix]) > _eps:
        if (k[ix] > 0):
            tx1 = (-0.5*w - r[:, ix])/k[ix]
            tx2 = (0.5*w - r[:, ix])/k[ix]
        else:
            tx2 = (-0.5*w - r[:, ix])/k[ix]
            tx1 = (0.5*w - r[:, ix])/k[ix]
        res = [tx1, tx2]
    else:
        n = r.shape[0]
        uno = np.ones((n, 1))
        res = [-_inf*uno, _inf*uno]
    return res


def crossRadial(rad, rsq, ksq, rk):
    """Calculate cross-times in radial coordinates

    Arguments:
        rad -- radius
        rsq -- array of r^2 values (r = position)
        ksq -- k^2 (k = direction vector)
        rk -- array of r.dot(k) values
    Returns:
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
    """
    D2 = rk*rk - ksq*(rsq-rad**2)
    mask = np.array(D2 > 0, dtype=int)
    D = np.sqrt(np.abs(D2))
    tin = (-rk - D)/ksq*mask
    tout = (-rk + D)/ksq*mask
    return [tin, tout]


def crossCylinder(rad, r, k):
    """Calculate cross-times with cyllindric surface (axis along r[:,1] )

    Arguments:
        rad -- radius
        r -- array[:,3] of position vectors
        k -- direction vector
    Returns:
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
    """
    r1 = r[:, 0::2]
    k1 = k[0::2]
    rsq = np.sum(r1**2, axis=1)
    ksq = k1.dot(k1)
    rk = r1.dot(k1)
    res = crossRadial(rad, rsq, ksq, rk)
    return res


def crossSphere(rad, r, k):
    """Calculate cross-times with spherical surface

    Arguments:
        rad -- radius
        r -- array[:,3] of position vectors
        k -- direction vector
    Returns:
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
    """
    rsq = np.sum(r**2, axis=1)
    ksq = k.dot(k)
    rk = r.dot(k)
    res = crossRadial(rad, rsq, ksq, rk)
    return res


def crossShell(R1, R2, r, k):
    """Calculate cross-times with hollow sphere

    Arguments:
        R1, R2 -- radii, R1 < R2
        r -- position vectors
        k -- direction vector
    Returns:
        list of [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
        tin, tout are sorted on axis=1
    """
    [ti1, tf1] = crossSphere(R1, r, k)
    [ti2, tf2] = crossSphere(R2, r, k)
    in1 = np.array(tf1 > ti1, dtype=int)
    in2 = np.array(tf2 > ti2, dtype=int)
    inf1 = (1-in1)*_inf
    inf2 = (1-in2)*_inf
    tin = np.array([ti2*in2 + inf2, tf1*in1 + inf1]).T
    tout = np.array([ti1*in1 + inf1, tf2*in2 + inf2]).T
    tin.sort(axis=1)
    tout.sort(axis=1)
    return [tin, tout]


def crossShellCyl(R1, R2, r, k):
    """Calculate cross-times with hollow cylinder

    Arguments:
        R1, R2 -- radii, R1 < R2
        r -- position vectors
        k -- direction vector
    Returns:
        list of [tin, tout]
        where tin, tout are arrays of entry and exit times (tin <= tout)
        tin, tout are sorted on axis=1
    """
    [ti1, tf1] = crossCylinder(R1, r, k)
    [ti2, tf2] = crossCylinder(R2, r, k)
    in1 = np.array(tf1 > ti1, dtype=int)
    in2 = np.array(tf2 > ti2, dtype=int)
    inf1 = (1-in1)*_inf
    inf2 = (1-in2)*_inf
    tin = np.array([ti2*in2 + inf2, tf1*in1 + inf1]).T
    tout = np.array([ti1*in1 + inf1, tf2*in2 + inf2]).T
    tin.sort(axis=1)
    tout.sort(axis=1)
    return [tin, tout]


def quad2D(surf, z0, rho, size, r, k):
    """ Calculate cross-section with 2D curved surface

    The surface eq. is z = z0 + 0.5*(rho[0]*x**2 + rho[1]*y**2)
    Written in vector notation to speed up
    Arguments:
        surf -- surface sign
        rho -- array(2), curvatures along x and y
        size --  array(2), dimensions (x,y)
        r --  array(:,3), coordinates of the starting point
        k --  array(3), velocity vector
    Returns:
        [t1, sgn1] if there is only 1 solution
        [t1, sgn1, t2, sgn2] if there are 2 solutions
        or [] if there is no solution
        t1[:], t2[:] = cross section times for r[:]
        sgn[:] = sign(surf*N*k) for r[:]
              N = surface normal at the cross-section
              N is chosen with N[2] is positive
        sgn = 0 for rows of r with no solution
        Excludes solutions with (x,y) outside the size limits
        Excludes tangential solutions
    """
    n = r.shape[0]
    rr = r**2
    a = 0.5*(rho[0]*k[0]**2 + rho[1]*k[1]**2)
    b = rho[0]*k[0]*r[:, 0] + rho[1]*k[1]*r[:, 1] - k[2]
    c = 0.5*(rho[0]*rr[:, 0] + rho[1]*rr[:, 1]) - r[:, 2] + z0
    res = []

#    wt = size/k0
    if (abs(a) < _eps):  # plane
        uno = np.ones(n)
        n = b.size
        t = -c/b
        # rt = r + k*t
        N = np.array([0., 0., 1.])
        sgn = np.sign(surf*N.dot(k))
        res = [t, sgn*uno]
    else:
        uno = np.ones((n, 1))
        D2 = b**2 - 4.*a*c
        # apply condition D2>0
        iD2 = np.array((D2 > _eps), dtype=int)
        # this will exclude D2<=0 from solutions, but avoids NaN from sqrt
        D = np.sqrt(iD2*D2)
        rh = np.array([-rho[0], -rho[1], 1.])
        if (n > 0):
            # get times and positions of the cross points
            t1 = (-b - D)/(2*a)
            rt1 = r + k*t1.reshape(n, 1)
            # surface normals
            v = np.concatenate((rt1[:, 0:2], uno), axis=1)
            N = rh*v
            # signs of N * k, apply bundary and D2>0 conditions
            sgn1 = np.sign(surf*N.dot(k.T))*iD2
            # repeat for the 2nd solution
            t2 = (-b + D)/(2*a)
            rt2 = r + k*t2.reshape(n, 1)
            v = np.concatenate((rt2[:, 0:2], uno), axis=1)
            N = rh*v
            sgn2 = np.sign(surf*N.dot(k.T))*iD2
            res = [t1, sgn1, t2, sgn2]
    return res


def crossPaths(cm):
    """Calculate cross-sections through the shape iterior

    Arguments:
        cm -- array returned by the ShapeAbstract.cross method or its childs
    Returns:
        list [inside, ic, p1, p2]
        inside -- true if current position is inside the shape
        ic -- index of the layer containing or following the current position
          ("following" applies for a gap)
        p1 -- path preceding current position
        p2 -- path following current position
    """
    n = cm.shape
    if (cm.size == 0) or (n[1] < 2):
        return [False, 0, 0.0, 0.0]
    p1 = 0.
    p2 = 0.
    ic = -1
    inside = False
    # print('crossPaths ', cm)

    for i in range(n[0]):
        if (cm[i, 1] <= 0.0):
            p1 += cm[i, 1]-cm[i, 0]
        elif (cm[i, 0] >= 0.0):
            p2 += cm[i, 1]-cm[i, 0]
            if (ic < 0):
                ic = i
        else:
            p1 += -cm[i, 0]
            p2 += cm[i, 1]
            ic = i
            inside = True
    if (ic < 0):  # r must be behind the last layer
        ic = n[0]
    return [inside, ic, p1, p2]
