#=======================================================================
"""
model_llb.py
------------

Reduced-magnetization Landau–Lifshitz–Bloch (LLB) model
compatible with Two-Temperature Model (TTM).

Based on:
Raposo et al., Applied Sciences 10, 1307 (2020)

Author: Stan Daniels
Rewritten with numerical stability + physical consistency
"""
#=======================================================================
import numpy as np
from scipy.ndimage import laplace
#=======================================================================

class LandauLifshitzBlochModel2D:
    """
    2D reduced-magnetization LLB model.
    Magnetization m is dimensionless (|m| ~ 1).
    """

    #===================================================================
    def __init__(self, material, gamma=None):
        self.material = material
        self.Tc = material.Tc
        self.gamma = gamma if gamma is not None else material.gamma

        self.mu0 = 4*np.pi*1e-7  # H/m

    #===================================================================
    @staticmethod
    def m_equilibrium(T, Tc, beta=0.33):
        """
        Reduced equilibrium magnetization m_e(T).
        """
        T = np.asarray(T)
        return np.clip(1 - T / Tc, 0.0, 1.0)**beta
    
    #===================================================================
    @staticmethod
    def brillouin_prime(x, S=1/2):
        # For S = 1/2, B(x) = tanh(x)
        return 1.0 / np.cosh(x)**2
    
    #===================================================================
    def chi_parallel(self, Te):
        kB = 1.380649e-23
        muB = 9.274e-24

        m_e = self.m_equilibrium(Te, self.Tc)
        m_e = np.clip(m_e, 1e-6, None)

        x = (self.Tc / Te) * m_e
        Bp = self.brillouin_prime(x)

        chi = (
            self.mu0 * muB / (kB * Te)
            * Bp / np.maximum(1 - (self.Tc / Te) * Bp, 1e-6)
        )

        return chi

    #================= Field components ================================
    # ---------------- External Field ----------------
    def H_external(self, m, H_ext=None):
        Nx, Ny, _ = m.shape
        if H_ext is None:
            H_ext = np.zeros((Nx, Ny, 3))
        else:
            H_ext = np.broadcast_to(H_ext, (Nx, Ny, 3))
        return H_ext

    # ---------------- Exchange Field ----------------
    def H_exchange(self, m, dx=1.0):
        H_ex = np.stack([
            laplace(m[:,:,0], mode="reflect"),
            laplace(m[:,:,1], mode="reflect"),
            laplace(m[:,:,2], mode="reflect")
        ], axis=-1)
        H_ex *= 2 * self.material.A_ex / (self.mu0 * self.material.Ms0 * dx**2)
        return H_ex

    # ---------------- anisotropy Field ----------------
    def H_anisotropy(self, m, Te, n=2):
        ez = np.array([0.0, 0.0, 1.0])

        m_e = self.m_equilibrium(Te, self.Tc)
        m_e = np.clip(m_e, 1e-6, None)

        Ku_T = self.material.Ku * m_e**n

        H_an = (
            2 * Ku_T[..., None]
            / self.material.Ms0
        ) * m[..., 2, None] * ez * self.mu0

        return H_an


    # ---------------- Inverse Faraday Effect Field ----------------
    def H_IFE(self, m, laser_sequence=None, t=0.0, X=None, Y=None):
        H_IFE = np.zeros_like(m)
        if laser_sequence is not None and X is not None and Y is not None:
            H_IFE[...,2] = laser_sequence.total_IFE_field_grid(X, Y, t, self.material)
        return H_IFE

    # ---------------- Demagnetization Field ----------------
    def H_demag(self, m):
        """
        Compute the demagnetization field using a simple finite-difference approximation
        for a thin film (Nx, Ny, 3).
        Note: For more accuracy, FFT-based convolution with the demag tensor is recommended.
        """
        # Assuming uniform thickness along z and free boundary conditions
        m_xx = laplace(m[...,0], mode="reflect")
        m_yy = laplace(m[...,1], mode="reflect")
        m_zz = laplace(m[...,2], mode="reflect")

        # Simple isotropic approximation
        H_demag = - (m_xx[...,None] + m_yy[...,None] + m_zz[...,None]) * m
        return H_demag

    # ---------------- Longitudal Relaxation Field ----------------
    def H_longitudinal(self, m, Te):
        m_mag = np.linalg.norm(m, axis=-1)
        m_mag = np.clip(m_mag, 1e-12, None)

        m_e = self.m_equilibrium(Te, self.Tc)
        m_e = np.clip(m_e, 1e-6, None)

        chi = self.chi_parallel(Te)
        chi = np.clip(chi, 1e-6, None)

        Hm = np.zeros_like(m)

        below = Te < self.Tc
        above = ~below

        # T < Tc
        Hm[below] = (
            1.0 / (2 * chi[below])[..., None]
            * (1.0 - (m_mag[below]**2) / (m_e[below]**2))[..., None]
            * m[below]
        )

        # T > Tc
        Hm[above] = (
            -1.0 / chi[above][..., None]
            * (1.0 + (3*self.Tc / (5*(Te[above]-self.Tc))) * m_mag[above]**2)[..., None]
            * m[above]
        )

        return Hm

    #===================================================================
    def H_effective(self, m, Te, dx=1.0, H_ext=None, laser_sequence=None, t=0.0, X=None, Y=None):
        """
        Total effective field as sum of components.
        """
        H_ext_field = self.H_external(m, H_ext)
        H_ex_field = self.H_exchange(m, dx)
        H_an_field = self.H_anisotropy(m, Te)
        H_IFE_field = self.H_IFE(m, laser_sequence, t, X, Y)
        H_demag_field = self.H_demag(m)
        #H_long_field = self.H_longitudinal(m, Te)

        H_eff = H_ext_field + H_ex_field + H_an_field + H_IFE_field+ H_demag_field # + H_long_field 
        return H_eff
  
    #===================================================================
    def step_grid(
        self,
        m,
        Te,
        dt,
        dx=1.0,
        H_ext=None,
        laser_sequence=None,
        t=0.0,
        X=None,
        Y=None
    ):
        """
        Single explicit Euler LLB step.
        """

        H = self.H_effective(
            m, Te, dx=dx,
            H_ext=H_ext,
            laser_sequence=laser_sequence,
            t=t, X=X, Y=Y
        )

        alpha_par = np.clip(
            self.material.alpha_parallel_T(Te), 0.0, 1.0
        )[...,None]

        alpha_perp = np.clip(
            self.material.alpha_perp_T(Te), 0.0, 1.0
        )[...,None]

        gamma_eff = self.gamma / (1 + alpha_perp**2)

        m_norm = np.linalg.norm(m, axis=-1, keepdims=True)
        m_norm = np.clip(m_norm, 1e-12, None)

        # ---------------- LLB equation ----------------
        m_cross_H = np.cross(m, H)
        m_cross_m_cross_H = np.cross(m, m_cross_H)
        m_dot_H = np.sum(m * self.H_longitudinal(m, Te), axis=2)
        
        m_e = self.m_equilibrium(Te, self.Tc)[..., None]
        m_e = np.clip(m_e, 1e-6, None)

        dm_dt = (
            - gamma_eff * m_cross_H
            - gamma_eff * alpha_perp * m_cross_m_cross_H / m_norm**2
            + gamma_eff * alpha_par * m_dot_H[..., None] * m / m_norm**2
        )

        # ---------------- Euler update ----------------
        m_new = m + dt * dm_dt

        # ---------------- Safety limiter ----------------
        m_mag = np.linalg.norm(m_new, axis=-1, keepdims=True)
        m_new = np.where(
            m_mag > 2.0,
            2.0 * m_new / m_mag,
            m_new
        )

        return np.nan_to_num(m_new)
#=======================================================================
