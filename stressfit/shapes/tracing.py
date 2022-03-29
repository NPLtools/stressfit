# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:12:32 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
# from numpy.linalg import norm
_eps = 1.0e-12
_inf = 1.0e+12

def ell(theta, a, b, alpha):
    """Calculate x,y coordinates of a rotated ellipse.""" 
    co = np.cos(theta)
    si = np.sin(theta)
    x = a*co*np.cos(alpha) + b*si*np.sin(alpha)
    y = -a*co*np.sin(alpha) + b*si*np.cos(alpha)
    return [x,y]    

def ell_limits(a, b, alpha):
    """Calculate x,y limits of a rotated ellipse."""
    thx = np.arctan2(b*np.sin(alpha),a*np.cos(alpha))
    thy = np.arctan2(-b*np.cos(alpha),a*np.sin(alpha))
    limx = [ell(thx, a, b, alpha)[0],ell(thx+np.pi, a, b, alpha)[0]]
    limy = [ell(thy, a, b, alpha)[1],ell(thy+np.pi, a, b, alpha)[1]]
    return [max(limx), max(limy)]

def get_metric_ell2D(a, b, angle):
    """Calculate metric tensor for a 2D ellipse."""
    dia = np.diagflat([1/a**2,1/b**2])
    co = np.cos(angle)
    si = np.sin(angle)
    R = np.array([[co, si],[-si, co]])
    gm = R.dot(dia.dot(R.T))
    return gm


def crossPolygon(r, k, height, sides, centres, normals, axes):
    """Calculate cross-times with a polygon.
    
    Polygon walls are defined by the normal vectors aiming inside.
    No check is made if the walls form a complete object.
    
    Recipe:
    
        1. Calculate flight times to the walls.
        2. The signs of the time and r.dot(normal) decides if it is an 
           entry or exit.
        3. Entry is the maximum of entry times.
        4. Exit is the minimum of exit times.

    Parameters
    ----------
    r : array(:,3)
        position vectors
    k : array(3)
        velocity vectors
    height ; float
        Height of the object
    sides : array
        side lengths of teh walls
    centres : list of array[3]
        centres of the walls
    normals : list of array[3]
        normals to the walls

    Returns
    -------
    list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin < tout)
    """
    tin = -_inf*np.ones(r.shape[0])
    tout = _inf*np.ones(r.shape[0])
    for i in range(len(normals)):
        d = (r-centres[i]).dot(normals[i])
        vn = k.dot(normals[i])
        sgn = np.sign(vn)
        if abs(vn) > _eps:
            t = -d/vn
            kt = np.kron(t,k).reshape((-1,3))
            xw = (r + kt - centres[i]).dot(axes[i]) 
            yh = np.abs(r[:,1])
            ins = np.array((xw<0.5*sides[i]) & (yh<0.5*height), dtype=int)
            if sgn > 0:
                # directs inside the object
                tt = t*ins + _inf*(ins-1)
                tin = np.maximum(tin,tt)
            else:
                # directs outside the object
                tt = t*ins - _inf*(ins-1)
                tout = np.minimum(tout,tt)
            
    # inspect cross-section with top/bottom walls
    if abs(k[1]) > _eps:
        [ty1, ty2] = crossLayer(height, r, k, 1)
        tin = np.maximum(tin,ty1)
        tout = np.minimum(tout,ty2)
    return [tin, tout]
            
            
def crossLayer(w, r, k, ix):
    """Calculate cross-times with plane parallel layer.

    Parameters
    ----------
    w: float
        layer thickness
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vectors
    ix: int
        index of axis normal to the layer [0 .. 2]

    Returns
    -------
    list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin < tout)
    """
    if abs(k[ix]) > _eps:
        if (k[ix] > 0):
            tx1 = (-0.5*w - r[:, ix])/k[ix]
            tx2 = (0.5*w - r[:, ix])/k[ix]
        else:
            tx1 = (0.5*w - r[:, ix])/k[ix]
            tx2 = (-0.5*w - r[:, ix])/k[ix]
        res = [tx1, tx2]
    else:
        uno = np.ones((r.shape[0], 1))
        res = [-_inf*uno, _inf*uno]
    return res


def crossRadial(rad, rsq, ksq, rk):
    """Calculate cross-times with a circle or sphere.

    Parameters
    ----------
    rad: float
        radius
    rsq: array(:)
        r*r (sqrares of position coordinates)
    ksq: array(:)
        k*k (squares of velocity coordinates)
    rk: array(:)
        r.dot(k)
    
    Returns
    -------
    list [tin, tout]
        tin, tout are entry and exit times (tin < tout). 
        If tin>=tout, then there is no cross-section.
    """
    D2 = rk*rk - ksq*(rsq-rad**2)
    mask = np.array(D2 > 0, dtype=int)
    D = np.sqrt(np.abs(D2))
    tin = (-rk - D)/ksq*mask
    tout = (-rk + D)/ksq*mask
    return [tin, tout]



def crossCylinder(rad, r, k):
    """Calculate cross-times with cyllindric surface (axis along r[:,1]).

    Parameters
    ----------
    rad: float
        radius
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vectors
    Returns
    -------
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin = tout).
        If tin>=tout, then there is no cross-section.
    """
    r1 = r[:, 0::2]
    k1 = k[0::2]
    rsq = np.sum(r1**2, axis=1)
    ksq = k1.dot(k1)
    rk = r1.dot(k1)
    res = crossRadial(rad, rsq, ksq, rk)
    return res


def crossEllipse(gm, r, k):
    """Calculate cross-times with cyllindric surface with elliptic basis.
    
    Cylinder axis is along r[:,1].

    Parameters
    ----------
    gm: array
        metric tensor
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vectors
    Returns
    -------
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin = tout).
        If tin>=tout, then there is no cross-section.
    """
    r1 = r[:, 0::2]
    k1 = k[0::2]
    rsq = np.sum(r1*(r1.dot(gm)), axis=1)
    ksq = k1.dot(gm.dot(k1))
    rk = r1.dot(gm.dot(k1))
    res = crossRadial(1.0, rsq, ksq, rk)
    return res


def crossSphere(rad, r, k):
    """Calculate cross-times with spherical.

    Parameters
    ----------
    rad: float
        radius
    r: array(:,3)
        position vectors
    k: array(:,3)
        direction vectors
        
    Returns
    -------
        list [tin, tout]
        where tin, tout are arrays of entry and exit times (tin = tout).
        If tin>=tout, then there is no cross-section.
    """
    rsq = np.sum(r**2, axis=1)
    ksq = k.dot(k)
    rk = r.dot(k)
    res = crossRadial(rad, rsq, ksq, rk)
    return res


def crossShell(R1, R2, r, k):
    """Calculate cross-times with hollow sphere.

    Parameters
    ----------
    R1, R2: float
        inner and outer radii, R1<R2
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vector
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


def crossHollowCyl(R, ctr, r, k):
    """Calculate cross-times with hollow cylinder.

    Parameters
    ----------
    R: array_like(2)
        Inner and outer radii 
    ctr: array(2,3)
        coordinates of the centers of the inner and outer cylinder.
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vector
    Return
    ------
    list of [tin, tout]
        tin, tout are arrays of entry and exit times (tin <= tout).
        tin, tout are sorted on axis=1.
    """
    [ti1, tf1] = crossCylinder(R[0], r-ctr[0,:], k)
    [ti2, tf2] = crossCylinder(R[1], r-ctr[1,:], k)
    in1 = np.array(tf1 > ti1, dtype=int)
    in2 = np.array(tf2 > ti2, dtype=int)
    inf1 = (1-in1)*_inf
    inf2 = (1-in2)*_inf
    tin = np.array([ti2*in2 + inf2, tf1*in1 + inf1]).T
    tout = np.array([ti1*in1 + inf1, tf2*in2 + inf2]).T
    tin.sort(axis=1)
    tout.sort(axis=1)
    return [tin, tout]

def crossTubes(R, ctr, r, k):
    """Calculate cross-times with a cylinder with multiple coaxial holes.

    Parameters
    ----------
    R: array_like
        Inner and outer radii. R[0] is the outer radius, R[1:] are holes radii. 
    ctr: array(:,3)
        Coordinates of the centers of the inner and outer cylinders.
        The first elements is for the outer radius.
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vector
    Return
    ------
    list of [tin, tout]
        tin, tout are arrays of entry and exit times (tin <= tout).
        tin, tout are sorted on axis=1 (corresponds to positions r).
        Row (axis=0) correspond to the radii, R.  
    """
    nc = len(R)
    ti = nc*[0]
    tf = nc*[0]
    ins = nc*[0]
    inf = nc*[0]
    for i in range(nc):
        [ti[i], tf[i]] = crossCylinder(R[i], r-ctr[i,:], k)
        ins[i] = np.array(tf[i]>ti[i], dtype=int) 
        inf[i] = (1-ins[i])*_inf
    tin = nc*[0]
    tout = nc*[0]
    for i in range(1,nc):
        tin[i] = tf[i]*ins[i] + inf[i]
        tout[i] = ti[i]*ins[i] + inf[i]
    tin[0] = ti[0]*ins[0] + inf[0]
    tout[0] = tf[0]*ins[0] + inf[0]
    atin =  np.array(tin).T
    atout =  np.array(tout).T
    atin.sort(axis=1)
    atout.sort(axis=1)
    ins = np.array(atout>atin, dtype=int) 
    return [atin, atout, ins]

def crossETubes(gm, ctr, r, k):
    """Calculate cross-times with a cylinder with multiple coaxial holes.
    
    Like crossETubes, but assumes elliptic basis of the cylinder and holes.

    Parameters
    ----------
    gm: list of arrays
        Metric matrices. gm[0] is the outer surface, g[1:] for the holes. 
    ctr: array(:,3)
        Coordinates of the centers of the outer cylinder and inner inner holes.
        The first row is for the outer cylinder.
    r: array(:,3)
        position vectors
    k: array(3)
        velocity vector
    Return
    ------
    list of [tin, tout]
        tin, tout are arrays of entry and exit times (tin <= tout).
        tin, tout are sorted on axis=1 (corresponds to positions r).
        Row (axis=0) correspond to the radii, R.  
    """
    nc = len(gm)
    ti = nc*[0]
    tf = nc*[0]
    ins = nc*[0]
    inf = nc*[0]
    for i in range(nc):
        [ti[i], tf[i]] = crossEllipse(gm[i], r-ctr[i,:], k)
        ins[i] = np.array(tf[i]>ti[i], dtype=int) 
        inf[i] = (1-ins[i])*_inf
    tin = nc*[0]
    tout = nc*[0]
    for i in range(1,nc):
        tin[i] = tf[i]*ins[i] + inf[i]
        tout[i] = ti[i]*ins[i] + inf[i]
    tin[0] = ti[0]*ins[0] + inf[0]
    tout[0] = tf[0]*ins[0] + inf[0]
    atin =  np.array(tin).T
    atout =  np.array(tout).T
    atin.sort(axis=1)
    atout.sort(axis=1)
    ins = np.array(atout>atin, dtype=int) 
    return [atin, atout, ins]

def crossShellCyl(R1, R2, r, k):
    """Calculate cross-times with hollow cylinder.
    
    Deprecated. For backward compatibility only.
    """
    return crossHollowCyl([R1, R2], np.zeros((2,3)), r, k)


def quad2D(surf, z0, rho, size, r, k):
    """ Calculate cross-section with 2D curved surface.

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


#%% test
