"""
Plots the stellar mass as a function of redshift.
"""


import unyt
import os

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

from loadStellarMassObservationalData import read_obs_data

stellar_density_output = unyt.msun / (unyt.Mpc ** 3)


try:
    plt.style.use("mnras_durham")
except:
    pass

simulations = {
    "anarchy-du": r"Anarchy-DU",
}


def load_data(simulation):
    # First we need to figure out exactly how many snapshots
    # there are.

    filename = f"{simulation}/eagle"

    # Look for the number of files in the directory.
    i = 0
    while True:
        if os.path.isfile("{}_{:04d}.hdf5".format(filename, i)):
            i += 1
        else:
            break

        if i > 10000:
            raise FileNotFoundError(
                "Could not find the snapshots in the directory"
            )


    # Now we can actually load the data and parse it.
    stellar_densities = []
    scale_factors = []

    for snapshot in range(i + 1):
        data = load(f"{filename}_{:04d}.hdf5")

        try:
            stellar_mass = data.stars.masses
            total_stellar_mass = stellar_mass.sum()
            volume = data.metadata.boxsize[0] * data.metadata.boxsize[1] * data.metadata.boxsize[2]

            stellar_density = (total_stellar_mass / volume).to(stellar_density_output)

        except AttributeError:
            stellar_density = 0.0 * stellar_density_output

        stellar_densities.append(stellar_density)
        scale_factors.append(data.metadata.scale_factor)

    return stellar_densities, scale_factors


# Load the data for real

simulation_data = {k: load_data(k) for k in simulations.keys()}

observational_data = read_obs_data()

# Plot the interesting quantities
fig, ax = plt.subplots()

ax.loglog()

# Simulation data plotting

simulation_lines = []
simulation_labels = []

for simulation in simulation_data.keys():
    scale_factor, redshift, sfr = simulation_data[simulation]
    name = simulations[simulation]

    # High z-order as we always want these to be on top of the observations
    simulation_lines.append(
        ax.plot(scale_factor, sfr.value, label=name, zorder=10000)[0]
    )
    simulation_labels.append(name)

# Observational data plotting

observation_lines = []
observation_labels = []

for index, observation in enumerate(observational_data):
    if observation.old_eagle_result:
        observation_lines.append(
            ax.plot(
                observation.scale_factor,
                observation.stellar_density,
                label=observation.description,
                color="aquamarine",
                zorder=-10000,
                linewidth=1,
                alpha=0.5
            )[0]
        )
    else:
        observation_lines.append(
            ax.errorbar(
                observation.scale_factor,
                observation.stellar_density,
                yerr=observation.density_error,
                xerr=observation.scale_factor_error,
                label=observation.description,
                linestyle="none",
                marker="o",
                elinewidth=0.5,
                markeredgecolor="none",
                markersize=2,
                zorder=index # Required to have line and blob at same zodrer
            )
        )

    observation_labels.append(observation.description)


ax.set_xlabel("Redshift $z$")
ax.set_ylabel(r"Stellar Density $\rho_*$ [M$_\odot$ Mpc$^{-3}$]")


# Observational data
redshift_ticks = np.array([0.0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0])
redshift_labels = [
    "$0$",
    "$0.2$",
    "$0.5$",
    "$1$",
    "$2$",
    "$3$",
    "$5$",
    "$10$",
    "$20$",
    "$50$",
    "$100$",
]
a_ticks = 1.0 / (redshift_ticks + 1.0)

ax.set_xticks(a_ticks)
ax.set_xticklabels(redshift_labels)
ax.tick_params(axis="x", which="minor", bottom=False)

ax.set_xlim(1.02, 0.07)
ax.set_ylim(3e5, 4e9)


simulation_legend = ax.legend(
    simulation_lines, simulation_labels, markerfirst=False, loc=1, fontsize=6
)
observation_legend = ax.legend(
    observation_lines, observation_labels, markerfirst=True, loc=3, fontsize=6
)
ax.add_artist(simulation_legend)

fig.tight_layout()

fig.savefig("Stellar_Mass_Comparison.pdf")
