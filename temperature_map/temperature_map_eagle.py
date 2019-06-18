"""
Creates a temperature map of the snapshot that is provided
and saves it as temperature_map.png.
"""

import sys
import h5py

from swiftsimio import load
from swiftsimio.visualisation import scatter, kernel_gamma

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm

from tqdm import tqdm

import numpy as np

resolution = 2048

snapshot = sys.argv[1]
n_files = int(sys.argv[2])

try:
    output_filename = sys.argv[3]
except IndexError:
    output_filename = "temperature_map_eagle.png"

def load_data(filename: str, n_files: int, to_read: str):
    """
    Loads the data from snapshots made of multiple files.
    """

    output = []

    for file in tqdm(range(n_files), desc=f"Reading {to_read}"):
        current_filename = f"{filename}.{file}.hdf5"

        with h5py.File(current_filename, "r") as handle:
            output.append(handle[to_read][...])

    return np.concatenate(output)


with h5py.File(f"{snapshot}.0.hdf5", "r") as handle:
    boxsize = handle["Header"].attrs["BoxSize"]

x, y, _ = load_data(snapshot, n_files, "PartType0/Coordinates").T / boxsize
hsml = load_data(snapshot, n_files, "PartType0/SmoothingLength") / (kernel_gamma* boxsize)
temp = load_data(snapshot, n_files, "PartType0/Temperature")

# We need to construct sum(T_j W_ij) / sum(W_ij)
norm_grid = scatter(x, y, np.ones_like(temp), hsml, resolution)
temp_grid = scatter(x, y, temp, hsml, resolution)

# Imsave does not take a norm
normalized_grid = LogNorm(vmin=1e2, vmax=1e8)(temp_grid / norm_grid)

imsave(output_filename, normalized_grid, cmap="twilight")
