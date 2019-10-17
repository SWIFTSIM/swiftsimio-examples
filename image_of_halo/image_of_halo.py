"""
Creates images of a given (set) of halo(es).
"""

import h5py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from swiftsimio import load
from swiftsimio.visualisation.sphviewer import SPHViewerWrapper

from unyt import unyt_array, msun

halo_numbers_in_catalogue = [0]
catalogue_filename = "halo/halo.properties"
snapshot_filename = "eagle.hdf5"

# Distance to visualise at, relative to R_vir
distance_factor = 0.5
resolution = 1024


def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


images = {
    "dark_matter": {"masses": {"cmap": "inferno", "vmin": None, "vmax": None}},
    "gas": {
        "masses": {"cmap": "viridis", "vmin": None, "vmax": None},
        "temperatures": {"cmap": "twilight", "vmin": 1e4, "vmax": 1e8},
        "metal_mass_fractions": {
            "cmap": "cubehelix",
            "vmin": 1e-1 * 0.012,
            "vmax": 10 * 0.012,
        },
    },
    "stars": {"masses": {"cmap": "bone", "vmin": None, "vmax": None}},
}

# Do some pre-computation and set up the particles instances.
data = load(snapshot_filename)
run_name = data.metadata.header["RunName"].decode("utf-8")

instances = {
    ptype: {
        field: SPHViewerWrapper(getattr(data, ptype), smooth_over=field)
        for field in images[ptype].keys()
    }
    for ptype in images.keys()
}

# Need to load the catalogue data now, including halo positions and virial radii.
with h5py.File(catalogue_filename, "r") as handle:
    positions = unyt_array(
        [handle["Xc"][...], handle["Yc"][...], handle["Zc"][...]],
        units=data.units.length,
    ).T

    radii = unyt_array(handle["R_200mean"], units=data.units.length)

    masses = unyt_array(handle["Mass_200mean"], units=data.units.mass)

# We have now fully computed smoothing lengths, etc. for our runs.
for halo in halo_numbers_in_catalogue:
    x, y, z = positions[halo]
    radius = radii[halo] * distance_factor
    mass = f"${latex_float(masses[halo].to(msun).value)}$ M$_\odot$"

    # Now we can render the individual images out to file
    for ptype, field_dict in instances.values():
        for field, instance in field_dict.values():
            instance.get_camera(
                x=x,
                y=y,
                z=z,
                r=radius,
                t=0.0,
                p=0.0,
                zoom=1,
                roll=0.0,
                xsize=resolution,
                ysize=resolution,
                extent=None,
            )

            instance.get_scene()
            instance.get_render()

            fig, ax = plt.subplots(figsize=(8, 8))
            fig.subplots_adjust(0, 0, 1, 1)
            ax.axis("off")

            cmap = images[ptype][field]["cmap"]
            vmin = images[ptype][field]["vmin"]
            vmax = images[ptype][field]["vmax"]

            # Need to remove the mass weighting for properties other than surface density
            if field != "masses":
                mass_image = instances[ptype][field]["masses"].image
                unweighted = instance.image / mass_image
            else:
                unweighted = instance.image

            ax.imshow(
                unweighted.value,
                cmap=cmap,
                norm=LogNorm(vmin=vmin, vmax=vmax),
                origin="lower",
            )

            ax.text(
                0.975,
                0.975,
                f"$z={data.metadata.z:3.3f}$",
                color="white",
                ha="right",
                va="top",
                transform=ax.transAxes,
            )

            ax.text(
                0.025,
                0.975,
                f"Smoothed property: {ptype}.{field}",
                color="white",
                ha="left",
                va="top",
                transform=ax.transAxes,
            )

            ax.text(
                0.025,
                0.025,
                f"Halo mass: {mass}",
                color="white",
                ha="left",
                va="bottom",
                transform=ax.transAxes,
            )

            ax.text(
                0.975,
                0.025,
                f"Run: {run_name}",
                color="white",
                ha="right",
                va="bottom",
                transform=ax.transAxes,
            )

            fig.savefig(f"{halo}_{field}_{ptype}.png", dpi=resolution // 8)
            fig.close()

