"""
Plots the stellar birth densities.
"""

import unyt

import matplotlib.pyplot as plt
import numpy as np
import h5py

from tqdm import tqdm

birth_density_output_units = unyt.mh / (unyt.cm ** 3)
birth_density_output_units_print = "cm$^{-3}$"

try:
    plt.style.use("mnras_durham")
except:
    pass

simulations = {
    "/cosma5/data/Eagle/ScienceRuns/Planck1/L0012N0188/PE/EagleVariation_NOAGN/data/snapshot_015_z002p012/snap_015_z002p012": "EAGLE12 NoAGN"
}

number_of_files = 16
bins = np.logspace(-2, 4, 100) * unyt.mh / (unyt.cm ** 3)
bins.convert_to_units(birth_density_output_units)
bin_centers = [0.5 * (a + b) for a, b in zip(bins[:-1].value, bins[1:].value)]


def read_data(filename: str, n_files: int, to_read: str):
    """
    Loads the data from snapshots made of multiple files.
    """

    output = []

    for file in tqdm(range(n_files), desc=f"Reading {to_read}"):
        current_filename = f"{filename}.{file}.hdf5"

        with h5py.File(current_filename, "r") as handle:
            output.append(handle[to_read][...])

    return np.concatenate(output)


def load_data(simulation):
    birth_densities = read_data(simulation, number_of_files, "PartType4/BirthDensity")

    with h5py.File(f"{simulation}.0.hdf5", "r") as handle:
        cgs = handle["PartType4/BirthDensity"].attrs["CGSConversionFactor"]
        z = handle["Header"].attrs["Redshift"]
        h = handle["Header"].attrs["HubbleParam"]

    birth_densities *= cgs * unyt.g / (unyt.cm ** 3)

    birth_densities.convert_to_units(birth_density_output_units)

    hist, _ = np.histogram(birth_densities.value, bins.value, density=True)

    return hist, z


fig, ax = plt.subplots()

ax.loglog()

# Simulation data plotting

simulation_lines = []
simulation_labels = []

for simulation in simulations.keys():
    try:
        birth_density_probability_density, redshift = load_data(simulation)
        name = simulations[simulation]

        # High z-order as we always want these to be on top of the observations
        simulation_lines.append(
            ax.plot(
                bin_centers,
                birth_density_probability_density,
                label=name,
                zorder=10000,
                lw=1,
            )[0]
        )
        simulation_labels.append(name)
    except (FileNotFoundError, OSError):
        pass


ax.loglog()

ax.text(
    0.025, 0.975, f"$z={redshift:2.2f}$", transform=ax.transAxes, ha="left", va="top"
)

ax.set_xlabel(f"Stellar birth density [{birth_density_output_units_print}]")
ax.set_ylabel(r"Probability density")

ax.set_xlim(bin_centers[0], bin_centers[-1])
ax.set_ylim(1e-5, 1.0)


simulation_legend = ax.legend(
    simulation_lines, simulation_labels, markerfirst=False, loc=1, fontsize=4
)


fig.tight_layout()

fig.savefig("Eagle_stellar_birth_densities.pdf")
