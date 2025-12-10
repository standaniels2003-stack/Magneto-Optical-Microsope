#=======================================================================
"""
model_llb.py
------------


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
from scipy.ndimage import laplace
from .helper import Brillouin, Brillouin_prime
#=======================================================================
class LandauLifshitzBlochModel2D:
    """
    2D Landau-Lifshitz-Bloch (LLB) model compatible with TTM.
    """
    #===================================================================
    def __init__(self, material, gamma=None):
        self.material = material
        self.Ms0 = material.Ms0
        self.Tc = material.Tc
        self.gamma = gamma if gamma is not None else material.gamma

        # Physical constants
        self.mu0 = 4*np.pi*1e-7       # Vacuum permeability (H/m)
        self.muB = 9.274e-24          # Bohr magneton (J/T)
        self.kB = 1.381e-23           # Boltzmann constant (J/K)

    #===================================================================
    @staticmethod
    def Me_equilibrium(T, Tc, Ms0, beta=0.33):
        """
        Temperature-dependent equilibrium magnetization.
        Vectorized for arrays.
        """
        T = np.asarray(T)
        return Ms0 * np.clip(1 - T / Tc, 0, 1)**beta

    #===================================================================
    def H_effective(self, m_grid, Te_map, H_ext=None, dx=1.0, laser_sequence=None, t=0.0, X=None, Y=None):
        """
        Compute the total effective field H_eff (vectorized) for the entire grid.
        """
        Nx, Ny, _ = m_grid.shape

        # External field
        H_ext_array = np.zeros((Nx, Ny, 3)) if H_ext is None else np.broadcast_to(H_ext, (Nx, Ny, 3))

        # --- Exchange Field ---
        H_exch = np.stack([
            laplace(m_grid[:,:,0], mode='reflect') / dx**2,
            laplace(m_grid[:,:,1], mode='reflect') / dx**2,
            laplace(m_grid[:,:,2], mode='reflect') / dx**2
        ], axis=-1)
        H_exch *= 2 * self.material.A_ex / (self.mu0 * self.Ms0**2)

        # --- Anisotropy Field (uniaxial along z) ---
        ez = np.array([0,0,1])
        m_dot_ez = m_grid[..., 2]
        H_ani = (2 * self.material.Ku / self.Ms0**2) * m_dot_ez[..., None] * ez

        # --- Longitudinal relaxation field ---
        Me = self.Me_equilibrium(Te_map, self.Tc, self.Ms0)[..., None]
        x = self.Tc * Me[..., 0] / Te_map
        m_mag = np.linalg.norm(m_grid, axis=-1)
        m_mag = np.clip(m_mag, 1e-18, None)
        eps = 1e-18
        B_prime_safe = np.clip(Brillouin_prime(x), eps, None)

        chi_par =  (self.mu0 * self.muB / (self.kB * Te_map) * B_prime_safe) / (1 - self.Tc ** B_prime_safe / Te_map)

        factor = np.where(
            Te_map < self.Tc,
            1 / (2 * chi_par) * (1 - m_mag**2 / (Me[..., 0]**2)),
            -1 / chi_par * (1 + (3 * self.Tc * m_mag**2) / ( 5 * (Te_map - self.Tc)))
        )

        # Apply factor to magnetization vector
        H_m = factor[..., None] * m_grid
        

        # --- IFE Field ---
        H_IFE = np.zeros_like(m_grid)
        if laser_sequence is not None and X is not None and Y is not None:
            H_IFE[..., 2] = laser_sequence.total_IFE_field_grid(X, Y, t, self.material, m_grid)

        # --- Total H_eff ---
        H_eff = H_ext_array + H_exch + H_ani + H_IFE + H_m
        return H_eff

    #===================================================================
    def step_grid(self, m_grid, Te_map, dt, dx=1.0, H_ext=None, laser_sequence=None, t=None, X = None, Y = None):
        """
        Perform a single LLB time step on the entire 2D magnetization grid.
        Vectorized implementation.
        """
        H_eff = self.H_effective(m_grid, Te_map, H_ext=H_ext, dx=dx, laser_sequence=laser_sequence, t=t, X=X, Y=Y)
        alpha_par = np.clip(self.material.alpha_parallel_T(Te_map)[..., None], 0, 1)
        alpha_perp = np.clip(self.material.alpha_perp_T(Te_map)[..., None], 0, 1)
        gamma_prime = self.gamma / (1 + alpha_perp**2)

        m_norm = np.linalg.norm(m_grid, axis=-1) + 1e-18

        # LLB update (Euler step)
        m_cross_H = np.cross(m_grid, H_eff)
        m_cross_m_cross_H = np.cross(m_grid, m_cross_H)
        m_dot_H = np.sum(m_grid * H_eff, axis=-1)
        print(m_grid.shape)
        print(H_eff.shape)
        Me = self.Me_equilibrium(Te_map, self.Tc, self.Ms0, beta=0.36)[..., None]
        dm_dt = -gamma_prime * m_cross_H \
                - gamma_prime * alpha_perp * m_cross_m_cross_H / m_norm[:,:,None]**2 \
                + gamma_prime * alpha_par * m_dot_H / m_norm[:,:,None]**2

        m_new = m_grid + dm_dt * dt

        return np.nan_to_num(m_new)
#=======================================================================