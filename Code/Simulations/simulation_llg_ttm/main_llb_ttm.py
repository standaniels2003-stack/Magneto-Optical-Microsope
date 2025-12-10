#=======================================================================
"""
main_llb_ttm.py
-------------------
Ultrafast laser simulation on a ferromagnetic thin film.
- Integrates Two-Temperature Model (TTM) and Landau-Lifshitz-Bloch (LLB)
- Animates Te, Tl, laser power, H_eff magnitude, and magnetization (mz)
- Progress bars + minimal console logs for monitoring

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
from scipy.ndimage import gaussian_filter

from src.laser import LaserPulse, LaserSequence
from src.materials import FerromagneticMaterial
from src.model_ttm import TwoTemperatureModel2D
from src.model_llb import LandauLifshitzBlochModel2D

#=======================================================================
# --- Simulation start info ---
now = datetime.datetime.now()
print("==============================================")
print(f"Simulation started at: {now.strftime('%Y-%m-%d %H:%M:%S')}")
start_time = time.time()

#=======================================================================
# --- Spatial and temporal grids ---
Nx, Ny = 100, 100
x = np.linspace(-100e-6, 100e-6, Nx)
y = np.linspace(-100e-6, 100e-6, Ny)
XX, YY = np.meshgrid(x, y, indexing='ij')
t_grid = np.linspace(0, 300e-15, 300)
dt = t_grid[1] - t_grid[0]

#=======================================================================
# --- Material and laser setup ---
print("Initializing material and laser sequence...")
material = FerromagneticMaterial(
    name="FePt", Ms0=1.1e5, Tc=750, Ku=5e6, A_ex=1e-11,
    gamma=1.76e11, thickness=10e-9, Ce_coeff=1e3,
    Cl=3e6, ke=10, kl=50, gel=1e17, R=0.3, chi_IFE=5e-13
)

pulse_params = [
    (0e-6, 0e-6, -1, 100e-15),
    (0e-6, 0e-6, -1, 400e-15),
    (0e-6, 0e-6, -1, 800e-15),
    (0e-6, 0e-6, -1, 1200e-15),
    (0e-6, 0e-6, -1, 1600e-15),
    (0e-6, 0e-6, -1, 2000e-15),
    (0e-6, 0e-6, -1, 2400e-15),
    (0e-6, 0e-6, -1, 2800e-15),
    (0e-6, 0e-6, -1, 3200e-15),
    (0e-6, 0e-6, -1, 3600e-15),
    (0e-6, 0e-6, -1, 4000e-15),
]
pulses = [LaserPulse(F=2e-8, t_TM=material.thickness, tau_L=200e-15,
                      tau_delay=50e-15, x0=x0, y0=y0, d0=50e-6, sigma=sigma, t0=t0)
          for x0, y0, sigma, t0 in pulse_params]
sequence = LaserSequence(pulses)
print("Laser sequence initialized.")

#=======================================================================
# --- Initialize solvers ---
print("Initializing TTM and LLB solvers...")
ttm = TwoTemperatureModel2D(material, sequence, grid_shape=(Nx, Ny), dx=x[1]-x[0])
llb_solver = LandauLifshitzBlochModel2D(material)
m_grid = np.zeros((Nx, Ny, 3))
m_grid[:, :, 2] = material.Ms0
print("Initial magnetization set along +z.")

#=======================================================================
# --- Precompute max values for color scales ---
print("Precomputing laser power and H_eff max values...")
power_max = 0.0
Hmax = 0.0
for t_val in tqdm(t_grid, desc="Precomputing"):
    power_map = np.array([[sequence.total_power(xi, yj, t_val) for yj in y] for xi in x])
    power_max = max(power_max, np.max(power_map))

    # Simple test grid for H_eff (spins along +z)
    m_test = np.zeros((Nx, Ny, 3))
    m_test[:, :, 2] = 1.0
    H_eff_grid = llb_solver.H_effective(m_test, Te_map=np.full((Nx,Ny), 300), t=t_val, laser_sequence=sequence)
    Hmax = max(Hmax, np.max(np.linalg.norm(H_eff_grid, axis=-1)))

#=======================================================================
# --- Figure setup ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
(ax_power, ax_Heff, ax_mz), (ax_Te, ax_Tl, ax_empty) = axes
ax_empty.axis('off')

def setup_imshow(ax, title, vmin, vmax, cmap="inferno"):
    im = ax.imshow(np.zeros((Nx, Ny)), origin='lower', cmap=cmap,
                   extent=[x[0]*1e6, x[-1]*1e6, y[0]*1e6, y[-1]*1e6],
                   vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    fig.colorbar(im, ax=ax)
    return im

im_power = setup_imshow(ax_power, "Laser Power (W/m^2)", 0, power_max, cmap="hot")
im_Heff = setup_imshow(ax_Heff, "|H_eff| (T)", 0, Hmax, cmap="coolwarm")
im_mz = setup_imshow(ax_mz, "Magnetization mz", -material.Ms0, material.Ms0, cmap="bwr")
im_Te = setup_imshow(ax_Te, "Electron Temp (K)", 300, 800, cmap="inferno")
im_Tl = setup_imshow(ax_Tl, "Lattice Temp (K)", 300, 330, cmap="inferno")
plt.tight_layout()

#=======================================================================
# --- Animation update function ---
def update(frame):
    t_val = t_grid[frame]
    t_fs = t_val * 1e15 

    # Laser power map
    power_map = np.array([[sequence.total_power(xi, yj, t_val) for yj in y] for xi in x])

    # Step TTM
    ttm.step(power_map, dt)
    Te_map, Tl_map = ttm.get_temperatures()

    # Step LLB
    m_grid_new = llb_solver.step_grid(m_grid, Te_map, dt, dx=x[1]-x[0], laser_sequence=sequence, t=t_val, X=XX, Y=YY)
    np.copyto(m_grid, m_grid_new)

    # Smoothed H_eff for visualization
    H_eff_map = llb_solver.H_effective(m_grid, Te_map, laser_sequence=sequence, t=t_val, X=XX, Y=YY)
    H_eff_map_smooth = gaussian_filter(np.linalg.norm(H_eff_map, axis=-1), sigma=1)

    # Update images
    im_power.set_data(power_map.T)
    im_Heff.set_data(H_eff_map_smooth.T)
    im_mz.set_data(m_grid[:,:,2].T)
    im_Te.set_data(Te_map.T)
    im_Tl.set_data(Tl_map.T)

    # Update titles to show current time
    ax_power.set_title(f"Laser Power (W/m^2) at t = {t_fs:.0f} fs")
    ax_Heff.set_title(f"|H_eff| (T) at t = {t_fs:.0f} fs")
    ax_mz.set_title(f"Magnetization mz at t = {t_fs:.0f} fs")
    ax_Te.set_title(f"Electron Temperature (K) at t = {t_fs:.0f} fs")
    ax_Tl.set_title(f"Lattice Temperature (K) at t = {t_fs:.0f} fs")

    # Extract H_IFE by recomputing it separately
    H_IFE_map = np.zeros_like(m_grid)
    H_IFE_map[... ,2] = sequence.total_IFE_field_grid(XX, YY, t_val, material)

    # Periodic console log
    if frame % 50 == 0 or frame == len(t_grid)-1:
        tqdm.write(f"t = {t_val*1e15:.0f} fs, max mz = {np.max(m_grid[:,:,2]):.3f}, min mz = {np.min(m_grid[:,:,2]):.3f}")
        tqdm.write(f"t = {t_fs:.0f} fs, Te_max = {np.max(Te_map):.0f} K, H_IFE_min = {np.min(H_IFE_map):.2e} T")

    return im_power, im_Heff, im_mz, im_Te, im_Tl

#=======================================================================
# --- Save animation ---
script_dir = os.path.dirname(os.path.abspath(__file__))
animation_folder = os.path.join(script_dir, "animations")
os.makedirs(animation_folder, exist_ok=True)
output_path = os.path.join(animation_folder, "simulation_results_ttm_llb.mp4")

fps = 60
writer = FFMpegWriter(fps=fps, metadata=dict(artist='Stan Daniels'), bitrate=1800)
with writer.saving(fig, output_path, dpi=200):
    for i in tqdm(range(len(t_grid)), desc="Saving animation"):
        update(i)
        writer.grab_frame()

end_time = time.time()
h, rem = divmod(end_time - start_time, 3600)
m, s = divmod(rem, 60)
print(f"Animation saved to: {output_path}")
print(f"Total simulation time: {int(h)}h {int(m)}m {s:.2f}s")
