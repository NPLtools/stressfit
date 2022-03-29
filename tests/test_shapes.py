# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:53:09 2022

@author: saroun
"""
import numpy as np
import stressfit.shapes as shapes
import stressfit.graphs as gr


def test_polygon(verbose=False):
    _eps = 1e-10
    deg = np.pi/180
    sdir = [np.cos(60*deg), 0, np.sin(60*deg)]
    side = 10.0
    shape = shapes.create('PolygonBar', num=6, side=side, 
                          height=20, angle=30, sdir=sdir)    
    
    ki = np.array([0,0,1])
    kf = np.array([np.cos(30*deg), 0, np.sin(30*deg)])
    
    r = np.array([[0., 0., 0.],
                  [0., 0., -1.5*side],
                  [0.5*side, 0., -side],
                  [0.5*side, 0., 0],
                  [0.5*side, 0., side]
                  ])
    if verbose:
        gr.plotShape([30,30], 0, shape)
        gr.plotShape([30,30], 1, shape)
        gr.plotShape([30,30], 2, shape)
        
    
    # calculate key dimensions:
    a = 0.5*side*np.tan(30*deg)
    b = 0.5*side - a
    v = side*np.cos(30*deg) 
    
    # test cross
    #-------------    
    # predicted values
    tc = np.zeros((len(r),2))
    tc[0] = [-side, side]
    tc[1] = [0.5*side, 2.5*side]
    tc[2] = [a, a + 2*b + side]
    tc[3] = [-(b + 0.5*side), (b + 0.5*side)]
    tc[4] = [-side - 2*b - a, -a]        
    # test cross
    ans = shape.cross(r, ki)
    if verbose:
        print('Cross times')
        for i in range(len(ans)):
            print('r = [{:g},{:g},{:g}]\t: {:g}, {:g}'.format(*r[i], *ans[i]))    
    # evaluate
    for i in range(len(ans)):
        q = np.subtract(tc[i],ans[i])
        qry = np.all(abs(q)<=_eps)
        assert(qry)
    if verbose:
        print('cross passed')    
    
    # test depthLocal 
    #---------------
    # predicted values
    s = np.array(sdir)
    qd = np.zeros(len(r))
    for i in range(len(qd)):
        qd[i] = r[i].dot(s)     
    # test depthLocal
    [d, d2, ins] = shape.depthLocal(r)
    if verbose:
        print('Depth')
        for i in range(len(r)):
            print('r = [{:g},{:g},{:g}]\t: {:g}'.format(*r[i], d[i]))    
    # evaluate
    q = np.subtract(qd,d)
    qry = np.all(abs(q)<=_eps)
    assert(qry)
    if verbose:
        print('depthLocal passed')
    
    # test rayPaths 
    #---------------
    # predicted values
    p = np.zeros((len(r),2))
    p[0] = [side, side]
    p[1] = [0, 0]
    p[2] = [0, 0]
    p[3] = [b + 0.5*side, (v - 0.5*side)/np.cos(30*deg)]
    p[4] = [0, 0]
    p = np.array(p)
    # test rayPaths
    path = np.array(shape.rayPaths(r, ki, kf)).T
    if verbose:
        print('Paths')
        for i in range(len(r)):
            print('r = [{:g},{:g},{:g}]\t: {:g}, {:g}'.format(*r[i],*path[i]))
# evaluate
    for i in range(len(r)):
        q = ins[i]*np.subtract(p[i],path[i])
        qry = np.all(abs(q)<=_eps)
        assert(qry)
    if verbose:
        print('rayPaths passed')
    print('test_polygon OK')

#%%
test_polygon()

