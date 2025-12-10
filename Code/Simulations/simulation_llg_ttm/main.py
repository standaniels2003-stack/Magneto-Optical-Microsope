#=======================================================================
"""
main.py
-------
Ultrafast laser simulation on a ferromagnetic thin film.
- Uses LaserPulse, LaserSequence, and TwoTemperatureModel2D
- Computes laser power, IFE field, electron and lattice temperatures
- Plots temporal evolution and animations

Author: Stan Daniels
Date: 
"""
#=======================================================================
import os
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from matplotlib.animation import FFMpegWriter

from src.laser import LaserPulse, LaserSequence
from src.materials import FerromagneticMaterial
from src.model_ttm import TwoTemperatureModel2D
#=======================================================================
# --- Define constants---
mu0 = 4*np.pi*1e-7

#=======================================================================
# --- Print current time and date ---
now = datetime.datetime.now()

# Print nicely formatted
print("==============================================")
print(f"Simulation started at: {now.strftime('%Y-%m-%d %H:%M:%S')}")
start_time = time.time()

#=======================================================================
# --- Spatial and temporal grids ---
Nx, Ny = 100, 100
x = np.linspace(-100e-6, 100e-6, Nx)
y = np.linspace(-100e-6, 100e-6, Ny)
X, Y = np.meshgrid(x, y)

t_grid = np.linspace(0, 1800e-15, 1800)
dt = t_grid[1] - t_grid[0]

#=======================================================================
# --- Material and laser sequence ---
print("----------------------------------------------")
print("Initializing Material and laser sequence...")
material = FerromagneticMaterial(
    name="FePt", Ms0=1.1e6, Tc=750, Ku=5e6, A_ex=1e-11, gamma=1.76e11,
    thickness=10e-9, Ce_coeff=1e3, Cl=3e6, ke=10, kl=50,
    gel=1e17, R=0.3, chi_IFE=1e-23
)

pulses = []
pulse_params = [
    (0e-6, 0e-6, 1, 100e-15),
    (0e-6, 0e-6, -1, 200e-15),
    (50e-6, 50e-6, 1, 300e-15),
    (-50e-6, 50e-6, 1, 400e-15),
    (-50e-6, -50e-6, -1, 500e-15),
    (50e-6, -50e-6, -1, 600e-15)
]

for x0, y0, sigma, t0 in tqdm(pulse_params, desc="Creating laser pulses"):
    pulses.append(LaserPulse(F=2e-8, t_TM=material.thickness, tau_L=50e-15,
                             tau_delay=20e-15, x0=x0, y0=y0, d0=50e-6, sigma=sigma, t0=t0))

sequence = LaserSequence(pulses)
print("Completed")
#=======================================================================
# --- Initialize TTM solver ---
print("----------------------------------------------")
print("Initializing Two Temperature Model solver...")
# Wrap with tqdm just to show “progress” even if it's fast
for _ in tqdm(range(1), desc="TTM setup"):
    ttm = TwoTemperatureModel2D(material, sequence, grid_shape=(Nx, Ny), dx=x[1]-x[0])
print("Completed")
#=======================================================================
# --- Precompute color scales ---
print("----------------------------------------------")
print("Precomputing color scales...")
power_max = 0
IFE_max = 0

for xi in tqdm(x, desc="Power max computation"):
    for yi in y:
        for t in t_grid:
            val = sequence.total_power(xi, yi, t)
            if val > power_max:
                power_max = val
            val_IFE = sequence.total_IFE_field(xi, yi, t, material)
            if abs(val_IFE) > IFE_max:
                IFE_max = abs(val_IFE)

Te_min, Te_max = 300, 800
Tl_min, Tl_max = 300, 330
print("Completed")
#=======================================================================
# --- Set up figure and axes ---
print("----------------------------------------------")
print("Setting up figure and axes...")
for _ in tqdm(range(1), desc="Creating figure and axes"):
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    (ax_power, ax_IFE), (ax_Te, ax_Tl) = axes
print("Completed")
# --- Initial maps ---
print("----------------------------------------------")
print("Initializing map images...")
# Wrap initialization with tqdm for visual effect
for _ in tqdm(range(1), desc="Creating map images"):
    im_power = ax_power.imshow(np.zeros((Nx, Ny)), origin='lower', cmap='hot',
                               extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6],
                               vmin=0, vmax=power_max)
    ax_power.set_title("Laser Power (W/m^2)")
    ax_power.set_xlabel("x (um)")
    ax_power.set_ylabel("y (um)")
    cbar_power = fig.colorbar(im_power, ax=ax_power)
    cbar_power.set_label("W/m^2")

    im_IFE = ax_IFE.imshow(np.zeros((Nx, Ny)), origin='lower', cmap='coolwarm',
                           extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6],
                           vmin=-IFE_max, vmax=IFE_max)
    ax_IFE.set_title("IFE Field (T)")
    ax_IFE.set_xlabel("x (um)")
    ax_IFE.set_ylabel("y (um)")
    cbar_IFE = fig.colorbar(im_IFE, ax=ax_IFE)
    cbar_IFE.set_label("T")

    im_Te = ax_Te.imshow(np.zeros((Nx, Ny)), origin='lower', cmap='inferno',
                         extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6],
                         vmin=Te_min, vmax=Te_max)
    ax_Te.set_title("Electron Temperature (K)")
    ax_Te.set_xlabel("x (um)")
    ax_Te.set_ylabel("y (um)")
    cbar_Te = fig.colorbar(im_Te, ax=ax_Te)
    cbar_Te.set_label("K")

    im_Tl = ax_Tl.imshow(np.zeros((Nx, Ny)), origin='lower', cmap='inferno',
                         extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6],
                         vmin=Tl_min, vmax=Tl_max)
    ax_Tl.set_title("Lattice Temperature (K)")
    ax_Tl.set_xlabel("x (um)")
    ax_Tl.set_ylabel("y (um)")
    cbar_Tl = fig.colorbar(im_Tl, ax=ax_Tl)
    cbar_Tl.set_label("K")

plt.tight_layout()
print("Completed")
#=======================================================================
# --- Animation update function ---
def update(frame):
    t_val = t_grid[frame]
    
    # Compute instantaneous maps
    power_map = np.array([[sequence.total_power(xi, yi, t_val) for yi in y] for xi in x])
    IFE_map   = np.array([[sequence.total_IFE_field(xi, yi, t_val, material) for yi in y] for xi in x])
    
    # Step TTM
    ttm.step(power_map, dt)
    Te_map, Tl_map = ttm.get_temperatures()
    
    # Update images
    im_power.set_data(power_map.T)
    im_IFE.set_data(IFE_map.T)
    im_Te.set_data(Te_map.T)
    im_Tl.set_data(Tl_map.T)
    
    ax_power.set_title(f"Laser Power (W/m^2) at t={t_val*1e15:.0f} fs")
    ax_IFE.set_title(f"IFE Field (T) at t={t_val*1e15:.0f} fs")
    ax_Te.set_title(f"Electron Temperature (K) at t={t_val*1e15:.0f} fs")
    ax_Tl.set_title(f"Lattice Temperature (K) at t={t_val*1e15:.0f} fs")
    
    return im_power, im_IFE, im_Te, im_Tl

#=======================================================================
# --- Save animation with progress bar ---
script_dir = os.path.dirname(os.path.abspath(__file__))
animation_folder = os.path.join(script_dir, "animations")
os.makedirs(animation_folder, exist_ok=True)
output_path = os.path.join(animation_folder, "simulation_results_ttm.mp4")

fps = 30
writer = FFMpegWriter(fps=fps, metadata=dict(artist='Stan Daniels'), bitrate=1800)
print("----------------------------------------------")
with writer.saving(fig, output_path, dpi=200):
    for i in tqdm(range(len(t_grid)), desc="Saving animation"):
        update(i)
        writer.grab_frame()

# Calculate total simulation time
end_time = time.time()
total_time = end_time - start_time
hours, rem = divmod(total_time, 3600)
minutes, seconds = divmod(rem, 60)

# Print final "We are done" message
print(f"Completed, animation saved to: {output_path}")
print(f"Total simulation time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")
print("----------------------------------------------")