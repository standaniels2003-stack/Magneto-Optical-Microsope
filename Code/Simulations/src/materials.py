#=======================================================================
"""
materials.py
------------
Material properties and physical constants for micromagnetic simulations
of ultrafast all-optical switching (AOS) in ferromagnetic thin films.

This module includes temperature-dependent parameters, thermal capacities,
magnetization saturation, and damping constants, based on experimental
and theoretical values for Co/Pt bilayers.

Equations implemented:
- Landau-Lifshitz-Bloch (LLB) model
- Two-temperature model (2TM) coupling parameters

Author: Stan Daniels
Date: 
"""
#=======================================================================
import numpy as np
#=======================================================================
# Physical Constants
#-----------------------------------------------------------------------
mu0 = 4 * np.pi * 1e-7          # Vacuum permeability [H/m]
kB = 1.380649e-23               # Boltzmann constant [J/K]
muB = 9.2740100783e-24          # Bohr magneton [J/T]
#=======================================================================
# Materials
#-----------------------------------------------------------------------
