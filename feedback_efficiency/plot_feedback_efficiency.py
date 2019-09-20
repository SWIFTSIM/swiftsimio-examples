"""
Plots the feedback efficiency for a given snapshot.
"""

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

try:
    plt.style.use("mnras_durham")
except:
    pass

simulations = {"/path/to/sim": "Simulation Name"}

snapshot = "eagle_0005.hdf5"

data = {k: load(f"{k}/{snapshot}") for k in simulations.keys()}

# Let's get the feedback efficiency max/min for the first
# one to get the bins we want to use.
first_data = data[list(simulations.keys())[0]]

feedback_min = float(
    first_data.metadata.parameters["EAGLEFeedback:SNII_energy_fraction_min"]
)
feedback_max = float(
    first_data.metadata.parameters["EAGLEFeedback:SNII_energy_fraction_max"]
)
feedback_bins = np.linspace(feedback_min, feedback_max, 50)
feedback_bin_centers = [
    0.5 * (x + y) for x, y in zip(feedback_bins[1:], feedback_bins[:-1])
]

current_redshift = first_data.metadata.redshift

# Can begin the plotting now

fig, ax = plt.subplots()

ax.semilogy()

# Simulation data plotting

simulation_lines = []
simulation_labels = []

for simulation in data.keys():
    this_data = data[simulation]
    these_feedback_fractions, _ = np.histogram(
        data.stars.feedback_energy_fractions, bins=feedback_bins
    )
    name = simulations[simulation]

    # High z-order as we always want these to be on top of the observations
    simulation_lines.append(
        ax.plot(feedback_bin_centers, these_feedback_fractions, label=name)[0]
    )
    simulation_labels.append(name)


ax.set_xlim(feedback_min, feedback_max)
ax.set_ylim(1, None)

ax.set_xlabel("Feedback energy fraction $f_E$")
ax.set_ylabel("Number of particles")

ax.legend()

fig.tight_layout()

fig.savefig("Compare_Feedback_EnergyFraction.pdf")
