"""
Creates a temperature map of the snapshot that is provided
and saves it as temperature_map.png.
"""

import sys
from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid

from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from p_tqdm import p_map

resolution = 1920

snapshot_start = int(sys.argv[1])
snapshot_end = int(sys.argv[2])
snapshot_name = sys.argv[3]

try:
    output_filename = sys.argv[4]
except IndexError:
    output_filename = "temperature_map.mp4"

snaps = list(range(snapshot_start, snapshot_end + 1))

data = [load(f"{snapshot_name}_{snapshot:04d}.hdf5") for snapshot in snaps]

def create_t_map(number):
    this_data = data[number]

    # We need to construct sum(T_j W_ij) / sum(W_ij)
    norm_grid = project_gas_pixel_grid(this_data, resolution, None)
    temp_grid = project_gas_pixel_grid(this_data, resolution, "temperatures")

    return temp_grid / norm_grid


maps = p_map(create_t_map, snaps)

fig, ax = plt.subplots(figsize=(1,1), dpi=resolution)
fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

image = ax.imshow(
    maps[0],
    norm=LogNorm(vmin=1e2, vmax=1e8),
    cmap="twilight"
)

def animate(number):
    image.set_array(maps[number])

    return image,


animation = FuncAnimation(fig, animate, snaps, interval=50)

animation.save("temperature_map.mp4", dpi=resolution)
