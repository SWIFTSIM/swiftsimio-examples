"""
Plots the star formation history.
"""

import unyt

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

from loadObservationalData import read_obs_data

sfr_output_units = unyt.msun / (unyt.year * unyt.Mpc ** 3)


try:
    plt.style.use("mnras_durham")
except:
    pass

simulations = {
    "anarchy-du": r"Anarchy-DU",
}


def load_data(simulation):
    filename = f"{simulation}/SFR.txt"
    snapshot = f"{simulation}/eagle_0000.hdf5"

    data = np.genfromtxt(filename).T

    initial_snapshot = load(snapshot)
    units = initial_snapshot.units
    boxsize = initial_snapshot.metadata.boxsize
    box_volume = boxsize[0] * boxsize[1] * boxsize[2]

    sfr_units = initial_snapshot.gas.star_formation_rates.units

    # a, Redshift, SFR
    return data[2], data[3], (data[7] * sfr_units / box_volume).to(sfr_output_units)


simulation_data = {k: load_data(k) for k in simulations.keys()}

observational_data = read_obs_data()

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
    if observation.fitting_formula:
        if observation.description == "EAGLE NoAGN":
            observation_lines.append(
                ax.plot(
                    observation.scale_factor,
                    observation.sfr,
                    label=observation.description,
                    color="aquamarine",
                    zorder=-10000,
                    linewidth=1,
                    alpha=0.5
                )[0]
            )
        else:
            observation_lines.append(
                ax.plot(
                    observation.scale_factor,
                    observation.sfr,
                    label=observation.description,
                    color="grey",
                    linewidth=1,
                    zorder=-1000,
                )[0]
            )
    else:
        observation_lines.append(
            ax.errorbar(
                observation.scale_factor,
                observation.sfr,
                observation.error,
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
ax.set_ylabel(r"SFR Density $\dot{\rho}_*$ [M$_\odot$ yr$^{-1}$ Mpc$^{-3}$]")


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
ax.set_ylim(1.8e-4, 1.7)


simulation_legend = ax.legend(
    simulation_lines, simulation_labels, markerfirst=False, loc=1, fontsize=6
)
observation_legend = ax.legend(
    observation_lines, observation_labels, markerfirst=True, loc=3, fontsize=6
)
ax.add_artist(simulation_legend)

fig.tight_layout()

fig.savefig("SFH_Comparison.pdf")
