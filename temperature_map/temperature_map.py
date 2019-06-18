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

try:
    output_filename = sys.argv[2]
except IndexError:
    output_filename = "temperature_map.png"

data = load(snapshot)

# We need to construct sum(T_j W_ij) / sum(W_ij)
norm_grid = project_gas_pixel_grid(data, resolution, None)
temp_grid = project_gas_pixel_grid(data, resolution, "temperature")

# Imsave does not take a norm
normalized_grid = LogNorm(vmin=1e2, vmax=1e8)(temp_grid / norm_grid)

imsave(output_filename, normalized_grid, cmap="twilight")
