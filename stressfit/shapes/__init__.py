"""Defines classes for various sample shapes.
    
All classes are dscendants of the abstract class shapeAbstract.


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
    `ShapeShellCyl(radius1, radius2, height)`
    
     - axis || y and outer/inner radii = radius1/radius2. 
          
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
from .shapeShellCyl import ShapeShellCyl
from .shapeSph import ShapeSph



