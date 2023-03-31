"""Implements classes for sample containers of various shapes. 
The class methods permit to calculate neutron ray paths, flight times to 
sample surfaces and depths under surface.
    
All classes are descendants of the abstract class shapeAbstract.

Use shapes.help() to print detailed information on individual shapes.

Examples
--------    
Infinite plate:
    `shapes.create(shapes.Plate,thickness=10)`

Cylinder:
    `shapes.create(shapes.Cylinder, radius=10)`

Sphere:
    `shapes.create(shapes.Sphere, radius=12)`

Tube:
    `shapes.create(shapes.Tube, Rin=4, Rout=12, height=30, ctr=[0,0], sref=0)`
     
Hollow sphere:
    `shapes.create(shapes.Shell, Rin=4, Rout=12)`    
    
Curved plate:
    `shapes.create(shapes.PlateCurved, Rin=4, Rout=12)`
    `shapes.create(shapes.PlateCurved, rho1=[0.03, 0.0], rho2=[0.03, 0.0])
    
    - z is along thickness, y is along height
    - rho1 = [hor, ver] are curvatures on the front surface
    - rho2 = [hor, ver] are curvatures on the rear surface
   
"""

from .shapeAbstract import ShapeAbstract
from .shapePlate import ShapePlate
from .shapeCyl import ShapeCyl
from .shapePlateCurved import ShapePlateCurved
from .shapeShell import ShapeShell
from .shapeTube import ShapeTube
from .shapeSph import ShapeSph
from .shapeTubes import ShapeTubes
from .shapeETubes import ShapeETubes
from .shapePolygonBar import ShapePolygonBar
import json

# constants for shape identification
# The strings should match the class names without the 'Shape' prefix
Plate = 'Plate'
PlateCurved = 'PlateCurved'
Cylinder = 'Cyl'
Tube = 'Tube'
Sphere = 'Sph'
Shell = 'Shell'
Tubes = 'Tubes'
ETubes = 'ETubes'
PolygonBar = 'PolygonBar'

def create(shape, **kwargs):
    """Create an instance of a shape class.
     
     Parameters
     ----------
     shape:
     One of the shape ID constant defined by this package:
     
     - stressfit.shapes.Plate 
     - stressfit.shapes.PlateCurved
     - stressfit.shapes.Cylinder
     - stressfit.shapes.Tube
     - stressfit.shapes.Sphere
     - stressfit.shapes.Shell
     - stressfit.shapes.Tubes
     - stressfit.shapes.ETubes
     - stressfit.shapes.PolygonBar
    
     **kwargs:
     Named arguments to the shape constructor.
     Use stressfit.shapes.help() to print detailed information.
     
    """
    comp = None
    try:
       if shape==Plate:
            comp = ShapePlate(**kwargs)
       elif shape==PlateCurved:
            comp = ShapePlateCurved(**kwargs)
       elif shape==Cylinder:
            comp = ShapeCyl(**kwargs)
       elif shape==Tube:
            comp = ShapeTube(**kwargs)
       elif shape==Sphere:
            comp = ShapeSph(**kwargs)
       elif shape==Shell:
            comp = ShapeShell(**kwargs)
       elif shape==Tubes:
            comp = ShapeTubes(**kwargs)
       elif shape==ETubes:
            comp = ShapeETubes(**kwargs)
       elif shape==PolygonBar:
            comp = ShapePolygonBar(**kwargs)
       elif shape=='File':
            comp = from_file(**kwargs)
       else:
           raise Exception('Undefined shape: {}'.format(shape))
    except Exception as e:
        print('Error in {}\n'.format(shape))
        print(e)
    return comp

def from_file(filename=''):
    """Create a shape instance form a file in JSON format. 
    
    The file structure should correspond to the output of the 
    :meth:`save` method of given class.
    

    Parameters
    ----------
    filename : str
        Input file in JSON format.

    Returns
    -------
    Instance of a shape class (descendant of ShapeAbstract)

    """
    if not filename:
        return None
    with open(filename, 'r') as f:
        lines = f.readlines()
    inp = json.loads('\n'.join(lines))
    if 'stressfit' not in inp:
        raise Exception('Not a Stressfit file.')
    if 'shape' not in inp['stressfit']:
        raise Exception('No shape defined in the input file.')
    data = inp['stressfit']['shape']
    clsname = data['class']
    if clsname.startswith('Shape'):
        clsname = clsname[5:]
    param =  data['param']  
    obj = create(clsname,**param)     
    return obj
    
def help():
    print('Call create(shape, **kwargs) to create an instance of sample shape.')
    fmt = '\n'+50*'-'+'\n'+'stressfit.shapes.{}\n'+50*'-'+'\n'
    print(fmt.format('create'))
    print(create.__doc__)
    print(fmt.format('Plate'))
    print(ShapePlate.__doc__)
    print(fmt.format('PlateCurved'))
    print(ShapePlateCurved.__doc__)
    print(fmt.format('Cylinder'))
    print(ShapeCyl.__doc__)
    print(fmt.format('Tube'))
    print(ShapeTube.__doc__)
    print(fmt.format('Sphere'))
    print(ShapeSph.__doc__)
    print(fmt.format('Shell'))
    print(ShapeShell.__doc__)
    print(fmt.format('Tubes'))
    print(ShapeTubes.__doc__)
    print(fmt.format('ETubes'))
    print(ShapeETubes.__doc__)
    print(fmt.format('PolygonBar'))
    print(ShapePolygonBar.__doc__)    


    






