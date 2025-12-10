#=======================================================================
"""
model_ttm.py
------------
Two-Temperature Model (TTM) solver for ultrafast laser heating of
ferromagnetic thin films.

This module provides a fully 2D, object-oriented TTM solver that evolves
electron and lattice temperatures over time under laser irradiation. It
is designed to integrate seamlessly with laser pulse sequences and 
material definitions, making it reusable for advanced micromagnetic 
simulations, including studies of all-optical switching.

Substrate cooling is treated via a characteristic time τ_subs, which
represents the heat exchange between the thin film lattice and the substrate.

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
#=======================================================================
class TwoTemperatureModel2D:
    """
    Two-temperature model (TTM) solver in 2D.

    Simulates the coupled dynamics of electron (Te) and lattice (Tl)
    temperatures in a thin film under laser excitation, using material
    parameters and a given laser pulse sequence.

    Attributes
    ----------
    material : object
        Material object providing thermal properties and electron heat capacity Ce(Te).
    laser_sequence : object
        Sequence of laser pulses applied to the system.
    Nx, Ny : int
        Number of spatial grid points in x and y directions.
    dx : float
        Spatial grid spacing [m].
    Te : np.ndarray
        Electron temperature array [K].
    Tl : np.ndarray
        Lattice temperature array [K].
    G : float
        Electron-lattice coupling constant [W/m^3/K].
    ke : float
        Electron thermal conductivity [W/m/K].
    kl : float
        Lattice thermal conductivity [W/m/K].
    Cl : float
        Lattice heat capacity [J/m^3/K].
    R : float
        Surface reflectivity (dimensionless, 0-1).
    thickness : float
        Film thickness [m].
    tau_subs : float
        Characteristic lattice-substrate cooling time [s].
    """

    #===================================================================
    def __init__(
        self, 
        material, 
        laser_sequence, 
        grid_shape: tuple[int, int], 
        dx: float, 
        initial_Te: float = 300.0, 
        initial_Tl: float = 300.0,
        tau_subs: float = 0.9e-9
    ):
        """
        Initialize the 2D Two-Temperature Model solver.

        Parameters
        ----------
        material : object
            Material object with attributes ke, kl, Cl, gel, R, thickness, 
            and a callable Ce(Te) returning electron heat capacity.
        laser_sequence : object
            LaserSequence object defining the applied laser pulses.
        grid_shape : tuple of int
            Number of points in (Nx, Ny) spatial grid.
        dx : float
            Spatial grid spacing [m].
        initial_Te : float, optional
            Initial electron temperature [K], default is 300 K.
        initial_Tl : float, optional
            Initial lattice temperature [K], default is 300 K.
        tau_subs : float, optional
            Characteristic lattice-substrate cooling time [s], default 0.9 ns.

        Raises
        ------
        ValueError
            If any input parameters are invalid (negative temperatures, dx <= 0, etc.).
        """
        # --- Sanity checks ---
        if not isinstance(grid_shape, (tuple, list)) or len(grid_shape) != 2:
            raise ValueError("grid_shape must be a tuple of two integers (Nx, Ny).")
        if dx <= 0:
            raise ValueError("dx must be positive.")
        if initial_Te < 0 or initial_Tl < 0:
            raise ValueError("Initial temperatures must be non-negative.")
        if tau_subs <= 0:
            raise ValueError("tau_subs must be positive.")
        if material.thickness <= 0:
            raise ValueError("Material thickness must be positive.")
        if not (0 <= material.R <= 1):
            raise ValueError("Reflectivity R must be between 0 and 1.")

        self.material = material
        self.laser_sequence = laser_sequence
        self.Nx, self.Ny = grid_shape
        self.dx = dx
        self.tau_subs = tau_subs  # substrate cooling time

        # Initialize temperature fields
        self.Te = np.full((self.Nx, self.Ny), initial_Te, dtype=float)
        self.Tl = np.full((self.Nx, self.Ny), initial_Tl, dtype=float)

        # Extract material properties
        self.G = material.gel
        self.ke = material.ke
        self.kl = material.kl
        self.Cl = material.Cl
        self.R = material.R
        self.thickness = material.thickness

    #===================================================================
    def step(self, power_map: np.ndarray, dt: float, T_sub: float = 300) -> None:
        """
        Advance the electron and lattice temperatures by one time step.

        Parameters
        ----------
        power_map : np.ndarray
            2D array of laser power absorbed at each grid point [W/m^2].
            Must have shape (Nx, Ny).
        dt : float
            Time step [s].
        T_sub : float, optional
            Substrate temperature [K], default is 300 K.

        Raises
        ------
        ValueError
            If power_map shape does not match grid or dt <= 0.
        RuntimeError
            If temperatures drop below absolute zero.
        """
        if dt <= 0:
            raise ValueError("Time step dt must be positive.")
        if power_map.shape != (self.Nx, self.Ny):
            raise ValueError(f"power_map shape {power_map.shape} does not match grid {(self.Nx, self.Ny)}.")

        # --- Laplacians ---
        lap_Te = laplace(self.Te, mode='reflect') / self.dx**2
        lap_Tl = laplace(self.Tl, mode='reflect') / self.dx**2

        # --- Electron heat capacity ---
        Ce = self.material.Ce(self.Te)
        if np.any(Ce <= 0):
            raise ValueError("Computed electron heat capacity Ce contains non-positive values.")

        # --- Absorption term ---
        absorption_term = (1 - self.R) * power_map / self.thickness

        # --- Update electron temperature ---
        dTe_dt = (-self.ke * lap_Te - self.G * (self.Te - self.Tl) + absorption_term) / Ce
        self.Te += dTe_dt * dt

        # --- Update lattice temperature with substrate cooling ---
        dTl_dt = (-self.kl * lap_Tl - self.G * (self.Tl - self.Te) - (self.Tl - T_sub)/self.tau_subs) / self.Cl
        self.Tl += dTl_dt * dt

        # --- Sanity check for temperatures ---
        if np.any(self.Te < 0) or np.any(self.Tl < 0):
            raise RuntimeError("Temperature field went below absolute zero.")

    #===================================================================
    def get_temperatures(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Retrieve the current electron and lattice temperature fields.

        Returns
        -------
        Te : np.ndarray
            Electron temperature array [K].
        Tl : np.ndarray
            Lattice temperature array [K].
        """
        return self.Te.copy(), self.Tl.copy()

    #===================================================================
    def total_energy(self) -> float:
        """
        Compute total energy density of the system [J/m^3].

        Returns
        -------
        float
            Sum of electron and lattice energies over the entire grid.
        """
        Ce = self.material.Ce(self.Te)
        return np.sum(Ce * self.Te + self.Cl * self.Tl)
