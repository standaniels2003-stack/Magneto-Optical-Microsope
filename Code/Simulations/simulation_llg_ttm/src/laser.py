#=======================================================================
"""
laser.py
--------
Object-oriented implementation of ultrafast laser pulses used in
micromagnetic simulations of all-optical switching (AOS) in
ferromagnetic thin films.

This module defines two primary classes:

1. LaserPulse:
    - Represents a single ultrafast laser pulse.
    - Implements Gaussian spatial and temporal profiles.
    - Computes maximum power, full space-time power, and the
      optically-induced magnetic field via the Inverse Faraday Effect (IFE).
    - Includes optional decay delay for the IFE field.

2. LaserSequence:
    - Handles sequences of laser pulses.
    - Computes the cumulative laser power and IFE field at any
      spatial-temporal point due to multiple pulses.

References
----------
Raposo, V., et al. (2020). Micromagnetic modeling of all-optical switching 
in ferromagnetic thin films: The role of inverse Faraday effect and 
magnetic circular dichroism. Applied Sciences, 10(4), 1307.

Author: Stan Daniels
Date: 
"""
#=======================================================================
import numpy as np
import warnings
from .helper import radial_distance
from .materials import FerromagneticMaterial
#=======================================================================
class LaserPulse:
    """
    Represents a single Gaussian laser pulse for micromagnetic simulations.

    Parameters
    ----------
    F : float
        Laser fluence in J/m^2.
    t_TM : float
        Thickness of the ferromagnetic sample (m).
    tau_L : float
        Full Width at Half Maximum (FWHM) of the pulse duration (s).
    tau_delay : float
        Optional delay for the decay of the IFE field (s).
    x0 : float
        x-coordinate of the laser spot center (m).
    y0 : float
        y-coordinate of the laser spot center (m).
    d0 : float
        Full Width at Half Maximum (FWHM) of the laser spot (m).
    sigma : int, optional
        Laser helicity for IFE calculation (+1 for right-handed, -1 for left-handed).

    Raises
    ------
    ValueError
        If fluence, thickness, pulse duration, or spot size are non-positive.
    TypeError
        If input types are invalid.

    Notes
    -----
    The laser power is calculated as:
        P0 = F / (t_TM * tau_L)
    The Gaussian spatial and temporal profiles are defined using
    standard FWHM-to-sigma conversions.

    The IFE field is computed assuming a helicity-dependent amplitude
    and an optional decay after the pulse peak.
    """

    def __init__(self, F, t_TM, tau_L, tau_delay, x0, y0, d0, sigma=1, t0=0.0):
        for name, val in zip(
            ["F", "t_TM", "tau_L", "tau_delay", "x0", "y0", "d0", "t0"],
            [F, t_TM, tau_L, tau_delay, x0, y0, d0, t0]
        ):
            if not isinstance(val, (int, float)):
                raise TypeError(f"{name} must be numeric (int or float).")
            
        if F <= 0:
            raise ValueError("Laser fluence F must be positive.")
        if t_TM <= 0:
            raise ValueError("Sample thickness t_TM must be positive.")
        if tau_L <= 0:
            raise ValueError("Pulse duration tau_L must be positive.")
        if d0 <= 0:
            raise ValueError("Laser spot size d0 must be positive.")
        if sigma not in (-1, 1):
            warnings.warn("Laser helicity sigma should be +1 or -1. Defaulting to +1.")
            sigma = 1
        if tau_delay < 0:
            warnings.warn("tau_delay is negative; this may produce unexpected IFE behavior.")

        self.F = F
        self.t_TM = t_TM
        self.tau_L = tau_L
        self.tau_delay = tau_delay
        self.x0 = x0
        self.y0 = y0
        self.d0 = d0
        self.sigma = sigma
        self.t0 = t0

    def P0(self):
        """
        Maximum instantaneous laser power (W/m^2) at the pulse peak.

        Returns
        -------
        float
            Maximum power of the laser pulse.
        """
        return self.F / (self.t_TM * self.tau_L)
    
    def spatial_profile(self, x, y):
        """
        Gaussian spatial profile of the laser pulse.

        Parameters
        ----------
        x, y : float
            Spatial coordinates (m) where the profile is evaluated.

        Returns
        -------
        float
            Relative spatial intensity at (x, y) (unitless, 0-1).
        """
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
            raise TypeError("x, y, and t must be numeric values.")
        
        r = radial_distance(x, y, self.x0, self.y0)
        return np.exp(-r**2 / (self.d0**2 / (4 * np.log(2))))
    
    def temporal_profile(self, t, t0):
        """
        Gaussian temporal profile of the laser pulse.

        Parameters
        ----------
        t : float
            Time at which the profile is evaluated (s).
        t0 : float
            Time of the pulse peak (s).

        Returns
        -------
        float
            Relative temporal intensity at time t (unitless, 0-1).
        """
        if not isinstance(t, (int, float)) or not isinstance(t0, (int, float)):
            raise TypeError("t and t0 must be numeric values.")
        
        return np.exp(-4*np.log(2) * (t - t0)**2 / self.tau_L**2)

    def temporal_profile_IFE(self, t, t0):
        """
        Temporal profile of the Inverse Faraday Effect (IFE) field,
        optionally including a decay delay after the pulse peak.

        Parameters
        ----------
        t : float
            Time at which the IFE field is evaluated (s).
        t0 : float
            Time of the pulse peak (s).

        Returns
        -------
        float
            Relative IFE field amplitude at time t (unitless, 0-1).
        """
        if not isinstance(t, (int, float)) or not isinstance(t0, (int, float)):
            raise TypeError("t and t0 must be numeric values.")
        
        tau_eff = self.tau_L if t < t0 else self.tau_L + self.tau_delay
        return np.exp(-4*np.log(2) * (t - t0)**2 / tau_eff**2)
    
    def power(self, x, y, t):
        """
        Full space-time laser power at a given point and time.

        Parameters
        ----------
        x, y : float
            Spatial coordinates (m).
        t : float
            Time (s).

        Returns
        -------
        float
            Instantaneous laser power at (x, y, t) in W/m^2.
        """
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)) or not isinstance(t, (int, float)):
            raise TypeError("x, y, and t must be numeric values.")
        
        return self.P0() * self.spatial_profile(x, y) * self.temporal_profile(t, self.t0)
    
    def IFE_field(self, x, y, t, material):
        """
        Computes the IFE field at (x, y, t) using the material's response.

        Parameters
        ----------
        x, y : float
            Spatial coordinates (m).
        t : float
            Time (s).
        material : FerromagneticMaterial
            Material object that defines the IFE susceptibility.

        Returns
        -------
        float
            Local IFE magnetic field (T).
        """
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)) or not isinstance(t, (int, float)):
            raise TypeError("x, y, and t must be numeric values.")
        if not isinstance(material, FerromagneticMaterial):
            raise TypeError("material must be an instance of FerromagneticMaterial")
        
        laser_power = self.power(x, y, t)
        envelope = self.temporal_profile_IFE(t, self.t0)
        return material.IFE_field(laser_power, m_local=1.0, sigma=self.sigma, temporal_envelope=envelope)


#=======================================================================
class LaserSequence:
    """
    Represents a sequence of Gaussian laser pulses.

    Parameters
    ----------
    pulses : list of LaserPulse, optional
        List of LaserPulse objects to initialize the sequence.

    Raises
    ------
    TypeError
        If any item in pulses is not a LaserPulse.
    
    Notes
    -----
    This class allows easy computation of the cumulative laser power
    and IFE field from multiple pulses at any spatial-temporal point.
    """

    def __init__(self, pulses=None):
        """
        Initialize the sequence. 
        pulses: list of LaserPulse objects
        """
        if not all(isinstance(p, LaserPulse) for p in pulses):
                raise TypeError("All elements in pulses must be instances of LaserPulse.")
        
        self.pulses = pulses if pulses is not None else []

    def add_pulse(self, pulse):
        """
        Adds a LaserPulse to the sequence.

        Parameters
        ----------
        pulse : LaserPulse
            The pulse to add.
        """
        if not isinstance(pulse, LaserPulse):
            raise TypeError("pulse must be an instance of LaserPulse.")
        self.pulses.append(pulse)

    def total_power(self, x, y, t):
        """
        Computes the total laser power at a point (x, y, t) from all pulses.

        Parameters
        ----------
        x, y : float
            Spatial coordinates (m).
        t : float
            Time (s).

        Returns
        -------
        float
            Total instantaneous laser power in W/m^2.
        """
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)) or not isinstance(t, (int, float)):
            raise TypeError("x, y, and t must be numeric values.")
        
        total = 0.0
        for pulse in self.pulses:
            total += pulse.power(x, y, t)
        return total
    

    def total_IFE_field(self, x, y, t, material):
        """
        Sum the IFE field of all pulses at (x, y, t) using the material response.

        Parameters
        ----------
        x, y : float
            Spatial coordinates (m).
        t : float
            Time (s).
        material : FerromagneticMaterial
            Material object that defines the IFE susceptibility.
        """
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)) or not isinstance(t, (int, float)):
            raise TypeError("x, y, and t must be numeric values.")
        if not isinstance(material, FerromagneticMaterial):
            raise TypeError("material must be an instance of FerromagneticMaterial")
    
        return sum(pulse.IFE_field(x, y, t, material) for pulse in self.pulses)
    
    def total_IFE_field_grid(self, X, Y, t, material, m_grid=None):
        """
        Vectorized calculation of the total IFE field over the entire grid.
        Includes Magnetic Circular Dichroism (MCD) → essential for real switching!
        Called automatically from model_llb.py
        """
        total = np.zeros_like(X)
        mz_local = m_grid[..., 2] if m_grid is not None else np.ones_like(X)

        for pulse in self.pulses:
            dx = X - pulse.x0
            dy = Y - pulse.y0
            r2 = dx*dx + dy*dy
            spatial = np.exp(-r2 / (pulse.d0**2 / (4 * np.log(2))))
            temporal = np.exp(-4*np.log(2) * (t - pulse.t0)**2 / pulse.tau_L**2)
            power = pulse.P0() * spatial * temporal

            # MCD effect: absorption depends on local mz and helicity
            absorption = material.absorption(mz_local, pulse.sigma)
            effective_power = power * absorption

            # IFE field contribution (sign from sigma)
            total += pulse.sigma * material.chi_IFE * effective_power

        return total
#=======================================================================