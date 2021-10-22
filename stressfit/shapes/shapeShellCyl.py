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

    def __init__(self, Rin=4.0, Rout=8.0, height=30.0):
        """Define hollow cylinder with axis || y."""
        super().__init__(Rin=Rin, Rout=Rout, height=height)
        warnings.warn("ShapeShellCyl is deprecated, use ShapeTube instead.", 
                      DeprecationWarning)
