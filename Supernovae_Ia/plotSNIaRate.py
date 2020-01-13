"""
Plots the star formation history.
"""
import matplotlib

matplotlib.use("Agg")

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

simulations = {"/cosma6/data/dp004/dc-nobe1/swift-colibre/examples/EAGLE_ICs/EAGLE_25_NOAGN_exp_loweff": r"Anarchy-DU"}


def load_data(simulation):
    filename = f"{simulation}/SNIa.txt"

    data = np.genfromtxt(filename).T
    
    default_SNIa_rate_conversion = 1.022690e-12

    a = (data[4] + data[5])/2.
    z = (data[6] + data[7])/2.

    # a, Redshift, SFR
    return a, z, data[11] * default_SNIa_rate_conversion 


simulation_data = {k: load_data(k) for k in simulations.keys()}

observational_data = read_obs_data()

fig, ax = plt.subplots()

ax.loglog()

# Simulation data plotting

simulation_lines = []
simulation_labels = []

for simulation in simulation_data.keys():
    scale_factor, redshift, SNIa_rate = simulation_data[simulation]
    name = simulations[simulation]

    # High z-order as we always want these to be on top of the observations
    simulation_lines.append(
        ax.plot(scale_factor, SNIa_rate, label=name, zorder=10000)[0]
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
                    observation.SNIa_rate,
                    label=observation.description,
                    color="aquamarine",
                    zorder=-10000,
                    linewidth=1,
                    alpha=0.5,
                )[0]
            )
        else:
            observation_lines.append(
                ax.plot(
                    observation.scale_factor,
                    observation.SNIa_rate,
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
                observation.SNIa_rate,
                observation.error,
                label=observation.description,
                linestyle="none",
                marker="o",
                elinewidth=0.5,
                markeredgecolor="none",
                markersize=2,
                zorder=index,  # Required to have line and blob at same zodrer
            )
        )
    observation_labels.append(observation.description)


ax.set_xlabel("Redshift $z$")
ax.set_ylabel(r"SNIa rate $[\rm yr^{-1} \cdot Mpc^{-3}]$")


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
ax.set_ylim(1e-5,2e-4)


simulation_legend = ax.legend(
    simulation_lines, simulation_labels, markerfirst=False, loc=1, fontsize=6
)
observation_legend = ax.legend(
    observation_lines, observation_labels, markerfirst=True, loc=3, fontsize=6, ncol=2
)
ax.add_artist(simulation_legend)

fig.tight_layout()

fig.savefig("SNIa_rate_Comparison.pdf")
