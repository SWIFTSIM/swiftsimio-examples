"""
Plots the feedback efficiency for a given snapshot, now using EAGLE.
"""

import matplotlib.pyplot as plt
import numpy as np

import read_eagle

try:
    plt.style.use("mnras_durham")
except:
    pass

# Point me towards ...../snapshot_x_z_x.0.hdf5
simulations = {"/path/to/sim": "Simulation Name"}


def load(snapshot_path):
    """
    Loads the data required.
    """

    snap = read_eagle.EagleSnapshot(snapshot_path)

    # Read whole volume...
    snap.select_region(0, 100, 0, 100, 0, 100)

    eff = snap.read_dataset(4, "Feedback_EnergyFraction")
    z = snap.redshift

    return z, eff


data = {k: load(k) for k in simulations.keys()}

# Let's get the feedback efficiency max/min for the first
# one to get the bins we want to use.
first_data = data[list(simulations.keys())[0]]

feedback_min = 0.3
feedback_bins = 3.0
feedback_bin_centers = [
    0.5 * (x + y) for x, y in zip(feedback_bins[1:], feedback_bins[:-1])
]

current_redshift = first_data[0]

# Can begin the plotting now

fig, ax = plt.subplots()

ax.semilogy()

# Simulation data plotting

simulation_lines = []
simulation_labels = []

for simulation in data.keys():
    this_data = data[simulation][1]
    these_feedback_fractions, _ = np.histogram(this_data, bins=feedback_bins)
    name = simulations[simulation]

    # High z-order as we always want these to be on top of the observations
    simulation_lines.append(
        ax.plot(feedback_bin_centers, these_feedback_fractions, label=name)[0]
    )
    simulation_labels.append(name)


ax.text(
    0.05,
    0.96,
    f"$z={current_redshift:2.2f}$",
    ha="left",
    va="top",
    transform=ax.transAxes,
)

ax.set_xlim(feedback_min, feedback_max)
ax.set_ylim(1, None)

ax.set_xlabel("Feedback energy fraction $f_E$")
ax.set_ylabel("Number of particles")

ax.legend()

fig.tight_layout()

fig.savefig("Compare_Feedback_EnergyFraction.pdf")
