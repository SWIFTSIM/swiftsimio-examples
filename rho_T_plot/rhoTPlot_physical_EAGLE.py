"""
Makes a rho-T plot. Uses the swiftsimio library.
"""
import matplotlib.pyplot as plt
import numpy as np 

import h5py

from unyt import mh, cm, Gyr
from tqdm import tqdm
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

# Constants; these could be put in the parameter file but are rarely changed.
density_bounds = [1e-7, 1e3]  # in nh/cm^3
temperature_bounds = [10**(3), 10**(8)]  # in K
bins = 128

# Plotting controls
cmap = "viridis"


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


def get_data(filename, n_files):
    """
    Grabs the data (T in Kelvin and density in mh / cm^3).
    """

    density = load_data(filename, n_files, "PartType0/Density")
    # Already in K
    temperature = load_data(filename, n_files, "PartType0/Temperature")

    with h5py.File(f"{filename}.0.hdf5", "r") as handle:
        print(handle["Header"].attrs["NumPart_Total"])
        scale_factor = handle["Header"].attrs["Time"]
        h = handle["Header"].attrs["HubbleParam"]
        unit_density = handle["Units"].attrs["UnitDensity_in_cgs"]

    mh = 1.6737236e-24
    unit_density /= mh

    print(len(density))

    density *= unit_density
    density *= h**2
    density /= scale_factor**3

    return density, temperature

def make_hist(filename, n_files, density_bounds, temperature_bounds, bins):
    """
    Makes the histogram for filename with bounds as lower, higher
    for the bins and "bins" the number of bins along each dimension.

    Also returns the edges for pcolormesh to use.
    """

    density_bins = np.logspace(
        np.log10(density_bounds[0]), np.log10(density_bounds[1]), bins
    )
    temperature_bins = np.logspace(
        np.log10(temperature_bounds[0]), np.log10(temperature_bounds[1]), bins
    )

    H, density_edges, temperature_edges = np.histogram2d(
        *get_data(filename, n_files), bins=[density_bins, temperature_bins]
    )

    return H.T, density_edges, temperature_edges


def setup_axes():
    """
    Creates the figure and axis object.
    """
    fig, ax = plt.subplots(1, figsize=(6, 5), dpi=300)

    ax.set_xlabel("Density [$n_H$ cm$^{-3}$]")
    ax.set_ylabel("Temperature [K]")

    ax.loglog()

    return fig, ax


def make_single_image(filename, n_files, density_bounds, temperature_bounds, bins):
    """
    Makes a single image and saves it to rhoTPlot_{filename}.png.
    
    Filename should be given _without_ hdf5 extension.
    """

    fig, ax = setup_axes()
    hist, d, T = make_hist(
        filename, n_files, density_bounds, temperature_bounds, bins
    )

    mappable = ax.pcolormesh(d, T, hist, cmap=cmap, norm=LogNorm())
    fig.colorbar(mappable, label="Number of particles", pad=0)

    fig.tight_layout()

    fig.savefig("rhoTPlot_EAGLE.png")

    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
             Plotting script for making a rho-T plot.
             Takes the filename handle, start, and (optionally) stop
             snapshots. If stop is not given, png plot is produced for
             that snapshot. If given, a movie is made.
             """
    )

    parser.add_argument(
        "-s",
        "--stub",
        help="""Stub for the filename (e.g. santabarbara). This is
                the first part of the filename for the snapshots,
                not including the final underscore. Required.""",
        type=str,
        required=True,
    )

    parser.add_argument(
         "-n",
         "--nfiles",
         help="""Number of EAGLE files to read""",
         type=int,
         required=True,
    )

    # Run in single image mode.
    args = vars(parser.parse_args())

    make_single_image(
        args["stub"],
        n_files=args["nfiles"],
        density_bounds=density_bounds,
        temperature_bounds=temperature_bounds,
        bins=bins,
    )
