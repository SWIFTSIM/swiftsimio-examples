"""
Creates a temperature map of the snapshot that is provided
and saves it as temperature_map.png.
"""

import sys
from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm

resolution = 2048

snapshot = sys.argv[1]

data = load(snapshot)

# We need to construct sum(T_j W_ij) / sum(W_ij)
norm_grid = project_gas_pixel_grid(snapshot, resolution, None)
temp_grid = project_gas_pixel_grid(snapshot, resolution, "temperature")

# Imsave does not take a norm
normalized_grid = LogNorm(vmin=1e2, vmax=1e8)(temp_grid / norm_grid)

imsave("temperature_map.png", normalized_grid, cmap="twilight")
