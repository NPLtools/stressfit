# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:53:09 2022

@author: saroun
"""
import numpy as np
import stressfit.shapes as shapes
import stressfit.dataio as dataio
import stressfit.graphs as gr



def _deep_comp(o1, o2):
    """Deep compare of two dictionaries or lists."""
    out = True
    try:
        if isinstance(o1, dict):
            k1 = o1.keys()
            k2 = o2.keys()
            assert len(k1)==len(k2)
            for key in k1:
                assert key in k2
            for key in k1:
                v1 = o1[key]
                v2 = o2[key]
                assert type(v1) == type(v2)
                if isinstance(v1, (dict, list)):
                    assert _deep_comp(v1, v2)
                else:
                    assert v1 == v2
        elif isinstance(o1, list):
            assert len(o1)==len(o2)
            for i in range(len(o1)):
                v1 = o1[i]
                v2 = o2[i]
                assert type(v1) == type(v2)
                if isinstance(v1, (dict, list)):
                    assert _deep_comp(v1, v2)
                else:
                    assert v1 == v2
        else:
            assert o1 == o2
    except:
        out = False
    return out
        
def test_shapes():     
    fn = dataio.get_output_file('tmp.json')    
#Infinite plate:
    obj = shapes.create(shapes.Plate,thickness=10)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Cylinder:
    obj = shapes.create(shapes.Cylinder, radius=10)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Sphere:
    obj = shapes.create(shapes.Sphere, radius=12)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Tube:
    obj = shapes.create(shapes.Tube, Rin=4, Rout=12, height=30, ctr=[0,0], 
                         sref=0)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Hollow sphere:
    obj = shapes.create(shapes.Shell, Rin=4, Rout=12)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Curved plate:
    obj = shapes.create(shapes.PlateCurved, rho1=[0.03, 0.0], 
                         rho2=[0.03, 0.0])
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#Tubes: 
    obj = shapes.create(shapes.Tubes, Rout=7, height=30, 
                         Rin=[3, 2], 
                         ctr=[[-4, 0],[3, 0]], 
                         sdir=[0,0,1],
                         sctr=[0,0,0])
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
#PolygonBar    
    obj = shapes.create('PolygonBar', num=6, side=8, 
                      height=20, angle=30)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    assert _deep_comp(obj.get_param(), obj2.get_param())
    print('stressfit.shapes.test OK')
    
def test_polygon(verbose=False, outfile='tmp.json'):
    _eps = 1e-10
    deg = np.pi/180
    sdir = [np.cos(60*deg), 0, np.sin(60*deg)]
    side = 10.0
    
    # define as regular polygon:
    #shape = shapes.create('PolygonBar', num=6, side=side, 
    #                      height=20, angle=30, sdir=sdir)  
    
    # define by edges
    edges = [
        [8.660254037844387, 5.0], 
        [8.660254037844389, -5.0], 
        [0.0, -10.0], 
        [-8.660254037844386, -5.0], 
        [-8.660254037844387, 5.0], 
        [0, 10.0]
        ]
    shape = shapes.create('PolygonBar', edges=edges, height=20, sdir=sdir) 
    fn = dataio.get_output_file(outfile)
    shape.save(fn)
    
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
    print('test PolygonBar OK')


def test_ETubes(outfile='tmp.json', verbose=False):
    def def_wires(param, inv=False):
        """Convert given list of wire parameters to a list of dict.
        
        Parameters
        ----------
        param : list
            A list to wire records. Each item is a list of numeric parameters for
            a single wire.
            - x, z : position of the centre
            - a, b : semi-axes of the elliptic basis
            - angle : angle of a with respect to x-axis
            
        Returns
        -------
        list
            The same list, but the items are converted to dict with 
            the parameter names.
        """
        keys = ['x', 'z', 'a', 'b', 'angle']
        wires = []
        for pset in param:
            h = {}
            for i in range(len(keys)):
                key = keys[i]
                h[key] =  pset[i]
            if inv:
                h['x'] = -h['x']
                h['angle'] = -h['angle']
            wires.append(h)
        return wires
    fn = dataio.get_output_file(outfile)
    # Define wire positions and shapes
    # dimensions are in [mm]
    # keys = ['x', 'z', 'a', 'b', 'angle']
    param = [[-0.59, -0.546, 3.6, 3.6, 0],
    [-0.329, -12.716, 3.25, 2.29, -5],
    [-0.339, 11.212, 3.41, 2.68, -17],
    [8.621, -10.672, 3.24, 2.17, -48],
    [-8.392, -8.105, 3.33, 2.41, 48]]
    
    wires = def_wires(param)
    obj = shapes.create(shapes.ETubes, a=17.5, b=17.5, height=50, 
                          holes=wires)
    obj.save(fn)
    obj2 = shapes.from_file(fn)
    wires2 = def_wires(param)
    obj2.update(holes=wires2)
    assert _deep_comp(obj.get_param(), obj2.get_param())
    if verbose:
        gr.plotShape([40,40], 1, obj)
    print('test ETubes OK')
    

#%%
test_shapes()
test_ETubes(verbose=False)
test_polygon(verbose=False)

"""

deg = np.pi/180
sdir = [np.cos(60*deg), 0, np.sin(60*deg)]
side = 10.0
shape = shapes.create('PolygonBar', num=6, side=side, 
                      height=20, angle=30, sdir=sdir)   
shape.save('shape.json')



#new_shape = shapes.from_file('shape.json')
#gr.plotShape([30,30], 1, new_shape)


new_shape2 = shapes.from_file('shape_edges.json')
gr.plotShape([30,30], 1, new_shape2)
new_shape2.save('shape_edges_output.json')
"""

