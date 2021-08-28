# -*- coding: utf-8 -*-
"""
Shape definition for a hollow cylinder
Created on Tue Aug 15 13:44:06 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""

from .shapeTube import ShapeTube
import warnings

class ShapeShellCyl(ShapeTube):
    """Define hollow cylinder with axis || y, coaxial.
    
    Equivalent to ShapeTube(Rin, Rout, height).
    For backward compatibility only. Use ShapeTube in new projects.

    """ 

    shape_type = 'cylinder_hollow'

    def __init__(self, Rin, Rout, height):
        """Define hollow cylinder with axis || y."""
        super().__init__(Rin, Rout, height)
        warnings.warn("ShapeShellCyl is deprecated, use ShapeTube instead.", 
                      DeprecationWarning)
