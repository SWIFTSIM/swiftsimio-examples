"""
Plots the stellar birth densities.
"""

import unyt

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

birth_density_output_units = unyt.mh / (unyt.cm ** 3)
birth_density_output_units_print = "cm$^{-3}$"

try:
    plt.style.use("mnras_durham")
except:
    pass

simulations = {
    k: k
    for k in [
        "AN_LIM_0p1_CONST_STD_STD",
        "AN_LIM_0p1_STD_NONE_STD",
        "AN_LIM_0p1_STD_STD_NOSMOOTH",
        "AN_LIM_0p1_STD_STD_STD",
        "AN_LIM_0p2_STD_STD_STD",
        "AN_NOLIM_0p1_STD_STD_STD",
        "AN_NOLIM_0p2_STD_STD_STD",
        # "MI_LIM_0p1_STD_STD_STD",
        # "MI_NOLIM_0p1_STD_STD_STD",
    ]
}

snapshot_number = 19
bins = np.logspace(-2, 4, 100) * unyt.mh / (unyt.cm ** 3)
bins.convert_to_units(birth_density_output_units)
bin_centers = [0.5 * (a + b) for a, b in zip(bins[:-1].value, bins[1:].value)]


def load_data(simulation):
    snapshot = f"{simulation}/eagle_{snapshot_number:04d}.hdf5"

    data = load(snapshot)

    data.stars.birth_densities.convert_to_units(birth_density_output_units)

    hist, _ = np.histogram(data.stars.birth_densities.value, bins.value, density=True)

    return hist, data.metadata.redshift


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

fig.savefig("Birth_density_comparison.pdf")
