"""
Creates a (scatter plot) video of a selected region, using
John Helly's EAGLE reading routines.

We will provide one of these for SWIFT soon.
"""

import yaml
import h5py
import read_eagle

import numpy as np

from p_tqdm import p_map

import matplotlib.pyplot as plt

try:
    plt.style.use("mnras_durham")
except:
    pass

from matplotlib.animation import FuncAnimation


def get_params(filename: str) -> dict:
    """
    Gets parameters from the yaml file and returns them as a dictionary.
    """

    with open(filename, "r") as handle:
        data = yaml.load(handle)

    return data


def periodic_in_more_than_one_dimension(region, boxsize) -> bool:
    """
    Checks if periodic in more than one dimension.

    region should have shape [[x0, x1], [y0, y1], [z0, z1]]
    """
    periodic = False

    for side in region:
        if side[0] < 0.0 or side[1] > boxsize:
            if periodic:
                return True
            else:
                periodic = True

    return False


def load_region(filename, region, particle_type, subsample, property="Coordinates"):
    """
    Loads a region (respecting periodic boundary conditions!)
    """

    with h5py.File(filename, "r") as handle:
        boxsize = handle["Header"].attrs["BoxSize"]
        redshift = handle["Header"].attrs["Redshift"]

    if periodic_in_more_than_one_dimension(region, boxsize):
        raise AttributeError(
            "Unable to process periodicity in more than one direction (load_region)"
        )

    snap = read_eagle.EagleSnapshot(filename)

    flat_region = np.array(region).flatten()

    try:
        # First load the particles that we 'know' to be in the region, ignoring
        # periodicity
        snap.select_region(*flat_region)
        read_data = snap.read_dataset(particle_type, property)

        # Now we need to deal with periodicity. First, find the periodic axis.

        for dimension, side in enumerate(region):
            if side[0] < 0.0:
                # Must be periodic on the 'left' side, here we wrap the box
                # around to the 'right' side
                new_region = flat_region
                new_region[dimension] = boxsize + new_region[dimension]
                new_region[dimension + 1] = boxsize
            elif side[1] > boxsize:
                # The periodic side is on the 'right' and we need to wrap back
                # around to the 'left'
                new_region = flat_region
                new_region[dimension] = 0.0
                new_region[dimension + 1] = side[1] - boxsize
            else:
                continue

            # Note that this can only, by definition above, happen once; no need to
            # worry about this happening multiple times and reading from file again
            snap.select_region(*new_region)
            data_new = snap.read_dataset(particle_type, property)

            # If our data is co-ordinates, we need to set the z-axis of that data
            # back to the original co-ordinate system.
            if property == "Coordinates":
                if side[0] < 0.0:
                    data_new[:, dimension] = data_new[:, dimension] - boxsize
                elif side[1] > boxsize:
                    data_new[:, dimension] = data_new[:, dimension] + boxsize

            # We _must_ break otherwise we'll keep adding or taking off the boxsize
            # (bad idea).

            read_data = np.concatenate([read_data, data_new])

            break


    except KeyError:
        read_data = None

    return read_data[::subsample], boxsize, redshift


def filename(snapshot_directory, snapshot_name):
    return f"{snapshot_directory}/{snapshot_name}/{snapshot_name.replace('shot', '')}.0.hdf5"


def read_data(params: dict):
    """
    Reads the data in parallel
    """

    filenames = [filename(params["snapshot_directory"], x) for x in params["snapshots"]]
    region = params["bounds"]

    def baked_reader(filename):
        # For passing to p_map
        return load_region(
            filename,
            region,
            params["particle_type"],
            params["subsample"],
            "Coordinates",
        )

    coordinate_data = p_map(baked_reader, filenames)

    coordinates = [x[0] for x in coordinate_data]
    redshifts = [x[2] for x in coordinate_data]

    try:

        def baked_reader(filename):
            # For passing to p_map
            return load_region(
                filename,
                region,
                params["particle_type"],
                params["subsample"],
                params["colour_by"],
            )[0]

        colour_data = p_map(baked_reader, filenames)
    except KeyError:
        colour_data = [np.ones(len(x)) for x in coordinates]

    return [c.T for c in coordinates], redshifts, colour_data


def make_video(params: dict) -> None:
    """
    Creates the video using the parameters in the dictionary and
    saves it to disk.
    """

    coordinates, redshifts, colours = read_data(params)
    region = params["bounds"]

    convert_direction = {"x": 0, "y": 1, "z": 2}
    horizontal = convert_direction[params["horizontal_axis"]]
    vertical = convert_direction[params["vertical_axis"]]

    fig, ax = plt.subplots(figsize=(4, 4))
    scatter_in, = ax.plot(
        coordinates[0][horizontal][colours[0] != 0],
        coordinates[0][vertical][colours[0] != 0],
        linestyle="none",
        ms=0.5,
        marker="o",
        markeredgecolor="none",
        color="C1",
    )

    scatter_out, = ax.plot(
        coordinates[0][horizontal][colours[0] == 0],
        coordinates[0][vertical][colours[0] == 0],
        linestyle="none",
        ms=0.5,
        marker="o",
        markeredgecolor="none",
        color="C2",
    )

    try:
        ax.axvline(params["vertical_line_at"], color="grey", ls="dashed", zorder=-1000, lw=1.0)
    except KeyError:
        pass

    try:
        ax.axhline(
            params["horizontal_line_at"], color="grey", ls="dashed", zorder=-1000, lw=1.0
        )
    except KeyError:
        pass

    ax.set_xlim(*region[horizontal])
    ax.set_ylim(*region[vertical])
    ax.set_xlabel(f"${params['horizontal_axis']}$ in simulation units [$h^{{-1}}$ Mpc]")
    ax.set_ylabel(f"${params['vertical_axis']}$ in simulation units [$h^{{-1}}$ Mpc]")

    ax.set_title(params["figure_title"])

    text = ax.text(
        0.95,
        0.95,
        f"$z={0.0:2.2f}$",
        va="top",
        ha="right",
        transform=ax.transAxes,
    )

    fig.tight_layout()

    def animate(snapshot):
        redshift = redshifts[snapshot]

        x_data = coordinates[snapshot][horizontal]
        y_data = coordinates[snapshot][vertical]
        colour_data = colours[snapshot]

        scatter_in.set_xdata(
                x_data[colour_data != 0]
        )
        scatter_in.set_ydata(
                y_data[colour_data != 0]
        )
        scatter_out.set_xdata(
                x_data[colour_data == 0]
        )
        scatter_out.set_ydata(
                y_data[colour_data == 0]
        )

        text.set_text(f"$z={redshift:2.2f}$")

        return (scatter_in, scatter_out)

    animation = FuncAnimation(
        # Your Matplotlib Figure object
        fig,
        # The function that does the updating of the Figure
        animate,
        # Frame information (here just frame number)
        np.arange(len(redshifts)),
        # Extra arguments to the animate function
        fargs=[],
        # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
        interval=1000 / 15,
    )

    animation.save(params["output_filename"])

    return


if __name__ == "__main__":
    from argparse import ArgumentParser as ap

    parser = ap(description="Creates a scatter plot movie")

    parser.add_argument(
        "-p", "--param", help="Parameter filename. Required", required=True, type=str
    )

    args = vars(parser.parse_args())

    params = get_params(args["param"])

    make_video(params=params)
