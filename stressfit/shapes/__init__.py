"""Implements classes for sample containers of various shapes. 
The class methods permit to calculate neutron ray paths, flight times to 
sample surfaces and depths under surface.
    
All classes are descendants of the abstract class shapeAbstract.


Examples
--------    
Infinite plate:
    `ShapePlate(thickness)`
    
    - z is along thickness

Cylinder:
    `ShapeCyl(radius, height)`
    
    - axis || y    
Sphere:
    `Sph(radius)`

Hollow cylinder:
    `ShapeTube([Rin, Rout], height, ctr=[0,0], sref=0)`
    
     - axis || y
     - ctr defines x,z position of the hole centre. 
     - sref defines the reference surface (from which depth is calculated),
       0 is for the inner surface
     
Hollow sphere:
    `ShapeShell(Rin, Rout)`      
    
Curved plate:
    `ShapePlateCurved(thickness, length, height, rho1, rho2)`
    
    - z is along thickness, y is along height
    - rho1 = [hor, ver] are curvatures on the front surface
    - rho2 = [hor, ver] are curvatures on the rear surface
   
"""


from .shapeAbstract import ShapeAbstract
from .shapeCyl import ShapeCyl
from .shapePlate import ShapePlate
from .shapePlateCurved import ShapePlateCurved
from .shapeShell import ShapeShell
from .shapeTube import ShapeTube
from .shapeSph import ShapeSph
from .shapeShellCyl import ShapeShellCyl


