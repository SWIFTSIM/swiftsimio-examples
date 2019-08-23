"""
Compares metallicilty in [N] runs, in both stars and gas.
Simply edit the 'simulations' and the 'file' that tells you
which snapshot and run to use.
"""

import numpy as np
import matplotlib.pyplot as plt
from swiftsimio import load

plt.style.use("mnras_durham")

file = "eagle_{SNAPSHOT_NUMBER}.hdf5"

simulations = {"/path/to/simulation": "Name On Plot"}

bins = np.logspace(-8, -1, 50)
plt.loglog()

for number, (path, name) in enumerate(simulations.items()):
    data = load(f"{path}/{file}")

    plt.hist(
        data.gas.metal_mass_fractions,
        bins=bins,
        label=f"Gas, {name}",
        histtype="step",
        color=f"C{number}",
    )
    plt.hist(
        data.stars.metal_mass_fractions,
        bins=bins,
        label=f"Stars, {name}",
        histtype="step",
        color=f"C{number}",
        linestyle="dotted",
    )

plt.xlabel("Metal Mass Fraction")

plt.text(
    0.05,
    0.95,
    f"$z={data.metadata.z:2.2f}$",
    transform=plt.gca().transAxes,
    ha="left",
    va="top",
)

plt.legend(fontsize=6, loc=(0.1, 0.4))
plt.tight_layout()

plt.savefig("metal_mass_fractions.pdf")
