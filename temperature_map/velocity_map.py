"""
Creates a temperature map of the snapshot that is provided
and saves it as velocity_map.png; this includes an overlay
of the velocity field.
"""

import sys
from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid, scatter

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm, Normalize

resolution = 2048

snapshot = sys.argv[1]

try:
    output_filename = sys.argv[2]
except IndexError:
    output_filename = "velocity_map.png"

data = load(snapshot)

# We need to construct sum(T_j W_ij) / sum(W_ij)
norm_grid = project_gas_pixel_grid(data, resolution, None)
temp_grid = project_gas_pixel_grid(data, resolution, "temperature")

# Imsave does not take a norm
normalized_grid = temp_grid / norm_grid

# Now we have to load the positions
boxsize = data.metadata.boxsize[0]
x, y, _ = data.gas.coordinates[:].T / boxsize
u, v, w = data.gas.velocities[:].T / data.gas.velocities.max()
h = data.gas.smoothing_length / boxsize

# Generate the smoothed velocity maps
u_grid = scatter(x, y, u, h, resolution).T / norm_grid.T
v_grid = scatter(x, y, v, h, resolution).T / norm_grid.T
w_grid = scatter(x, y, w, h, resolution).T / norm_grid.T
x = np.linspace(0, 1, resolution)
y = np.linspace(0, 1, resolution)

# Calculate the actual speed for the line width
speed = np.sqrt(u_grid * u_grid + v_grid * v_grid + w_grid * w_grid)
lw = 0.5 * speed / speed.max()

# Can't do imsave here, have to actually set up a figure.
fig, ax = plt.subplots(figsize=(8, 8), dpi=resolution // 8)
fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

ax.streamplot(
    x, y, u_grid, v_grid, density=10, arrowstyle="-", linewidth=lw, color="white"
)

ax.imshow(
    normalized_grid.T,
    cmap="twilight",
    norm=LogNorm(vmin=1e2, vmax=1e8),
    extent=[0, 1, 0, 1],
    origin="lower",
)

fig.savefig(output_filename)
