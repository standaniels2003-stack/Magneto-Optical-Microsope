#=======================================================================
"""
helper.py
---------
Utility functions for micromagnetic simulation: grid generation,
interpolation, spatial derivatives, and field calculations.

Author: Stan Daniels
Date: 
"""
#=======================================================================
import numpy as np
#=======================================================================
def radial_distance(x, y, x0=0.0, y0=0.0):
    """
    Compute the radial distance from a point (x0, y0).
    Works for scalars or numpy arrays.
    """
    return np.sqrt((x - x0)**2 + (y - y0)**2)