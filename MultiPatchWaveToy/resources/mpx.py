"""mpx: MultiPatchX plot and validation utility.

Usage:
  mpx augment-tsv <data-file> <coord-file>
  mpx plot-tsv <augmented-data> <var> [--rhs] [--save] [--diverging]
  mpx plot-openpmd <data-file> <thorn-name> <group-name> <var-name> <num-patches> <iteration> [--refinement-level=<level>] [--z-slice=<zval>] [--openpmd-format=<format>] [--save] [--diverging] [--verbose]
  mpx (-h | --help)
  mpx --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --rhs                       Plot RHS variables.
  --save                      Save plot to file.
  --verbose                   Show OpenPMD file information.
  --diverging                 Use a diverging color map.
  --refinement-level=<level>  The refinement level to plot [default: 0].
  --z-slice=<zval>            The (global) Z coordinate slice to plot [default: 0.0]
  --openpmd-format=<format>   The underlying format of OpenPMD files [default: .bp5].
"""

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import openpmd_api as io

import subprocess
import multiprocessing
import os

from docopt import docopt


def augment_tsv(state_file, coord_file):
    print("Augmenting", state_file, "with coordinate data from", coord_file)

    p1 = subprocess.Popen(
        [
            "cut",
            "-f12,13,14",
            coord_file
        ],
        stdout=subprocess.PIPE
    )

    fout = open(os.path.basename(state_file), "w")

    subprocess.run(
        [
            "paste",
            state_file,
            "-"
        ],
        stdin=p1.stdout,
        stdout=fout
    )

    p1.wait()
    fout.close()


def make_plot(x, y, z, var_name, save_file, diverging):
    plt.close("all")

    if diverging == True:
        plt.tricontourf(
            x,
            y,
            z,
            cmap="seismic",
            levels=100
        )
    else:
        plt.tricontourf(
            x,
            y,
            z,
            levels=100  # np.linspace(0.0, 0.04, 100),
        )

    cb = plt.colorbar()
    cb.ax.set_ylabel(var_name)

    plt.tight_layout()

    if not save_file:
        plt.show()
    else:
        plt.savefig(var_name + ".png")

    plt.close("all")


def plot_tsv(arguments):
    # Common file vars
    vars1 = [
        "iteration",
        "time",
        "patch",
        "level",
        "component",
        "i",
        "j",
        "k",
        "x",
        "y",
        "z",
    ]

    vars2 = [
        "vcoordx",
        "vcoordy",
        "vcoordz"
    ]

    # Arguments
    rhs = arguments["--rhs"]
    augmented_state_data = arguments["<augmented-data>"]
    var = arguments["<var>"]
    save_file = arguments["--save"]
    diverging = arguments["--diverging"]

    # Datafile vars, for pandas
    if rhs:
        vars = vars1 + ["u_rhs", "rho_rhs"] + vars2
    else:
        vars = vars1 + ["u", "rho"] + vars2

    print("Reading data")
    data = pd.read_csv(augmented_state_data,
                       delim_whitespace=True, names=vars, comment="#")

    # Z slice data
    print("Filtering data for the z = 0 slice")

    sliced_data = data.loc[
        (data["vcoordz"] == 0.0) &
        (np.abs(data["x"]) <= 1.0) &
        (np.abs(data["y"]) <= 1.0) &
        (np.abs(data["z"]) <= 1.0)
    ]

    # RHS error data
    if rhs:
        print("Creating error plot")

        expected = -3.0 * np.pi**2 * np.cos(np.pi * sliced_data["vcoordx"]) * np.cos(
            np.pi * sliced_data["vcoordy"])
        error = np.abs(expected - sliced_data[var])

        max_error = np.amax(error)
        error = error / max_error

        make_plot(sliced_data["vcoordx"], sliced_data["vcoordy"],
                  error, var + "_err", save_file, diverging)

    print("Creating data plot")
    make_plot(sliced_data["vcoordx"], sliced_data["vcoordy"],
              sliced_data[var], var, save_file, diverging)


def openpmd_to_dataframe(fname, verbose, iteration_index, mesh_name, gf_name):
    # Based on
    # https://gist.github.com/stevenrbrandt/ef3718d52680d8a4e6accbad198a0b42
    # with a few extra spices
    series = io.Series(fname, io.Access.read_only)

    print("Reading dataset:", fname)

    if iteration_index not in series.iterations:
        raise RuntimeError(
            "Iteration {} not available".format(iteration_index))

    if iteration_index not in series.iterations:
        raise RuntimeError(
            "Iteration {} not available".format(iteration_index))

    iter = series.iterations[iteration_index]
    uu = iter.meshes[mesh_name]
    data_raw = uu[gf_name]
    data = data_raw.load_chunk()

    series.flush()

    # Get Axis labels
    x_index = uu.axis_labels.index("x")
    y_index = uu.axis_labels.index("y")
    z_index = uu.axis_labels.index("z")

    # Get grid spacing
    dx = uu.grid_spacing[x_index]
    dy = uu.grid_spacing[y_index]
    dz = uu.grid_spacing[z_index]

    # Get Grid origin
    x0 = uu.grid_global_offset[x_index]
    y0 = uu.grid_global_offset[y_index]
    z0 = uu.grid_global_offset[z_index]

    chunk_dataframes = []

    for chunk in data_raw.available_chunks():
        chunk_i_0 = chunk.offset[x_index]
        chunk_j_0 = chunk.offset[y_index]
        chunk_k_0 = chunk.offset[z_index]

        chunk_nx = chunk.extent[x_index]
        chunk_ny = chunk.extent[y_index]
        chunk_nz = chunk.extent[z_index]

        chunck_i_range = range(chunk_i_0, chunk_i_0 + chunk_nx)
        chunck_j_range = range(chunk_j_0, chunk_j_0 + chunk_ny)
        chunck_k_range = range(chunk_k_0, chunk_k_0 + chunk_nz)

        # Rearrange data for dataframe
        chunk_data = []

        for k in chunck_k_range:
            for j in chunck_j_range:
                for i in chunck_i_range:
                    chunk_data.append([
                        i,
                        j,
                        k,
                        x0 + i * dx,
                        y0 + j * dy,
                        z0 + k * dz,
                        data[k, j, i]
                    ])

        chunk_df = pd.DataFrame(
            chunk_data,
            columns=["i", "j", "k", "x", "y", "z", gf_name]
        )
        chunk_dataframes.append(chunk_df)

    if verbose:
        print("Merging dataframes")

    merged_df = pd.concat(chunk_dataframes)

    return merged_df


def merge_multipatch_dataframes(verbose, fname, iteration_index, patch, level, mesh_name, gf_name, z_slice_value):
    if verbose:
        print("Loading coordinate data")

    vcoordx_df = openpmd_to_dataframe(
        fname,
        verbose,
        iteration_index,
        f"coordinatesx_vertex_coords_patch{patch}_lev{level}",
        "coordinatesx_vcoordx"
    )

    vcoordy_df = openpmd_to_dataframe(
        fname,
        verbose,
        iteration_index,
        f"coordinatesx_vertex_coords_patch{patch}_lev{level}",
        "coordinatesx_vcoordy"
    )

    vcoordz_df = openpmd_to_dataframe(
        fname,
        verbose,
        iteration_index,
        f"coordinatesx_vertex_coords_patch{patch}_lev{level}",
        "coordinatesx_vcoordz"
    )

    if verbose:
        print("Merging coordinate data")

    merge_cols = ["i", "j", "k", "x", "y", "z"]

    coords_df = pd.merge(vcoordx_df, vcoordy_df, on=merge_cols)
    coords_df = pd.merge(coords_df, vcoordz_df, on=merge_cols)

    if verbose:
        print("Loading main data")

    main_gf_df = openpmd_to_dataframe(
        fname,
        verbose,
        iteration_index,
        mesh_name,
        gf_name
    )

    if verbose:
        print("Merging main data with coordinates")
    merged_df = pd.merge(coords_df, main_gf_df, on=merge_cols)

    if verbose:
        print("Filtering main data for z = ", z_slice_value)
    filtered_df = merged_df[merged_df["coordinatesx_vcoordz"] == z_slice_value]

    return filtered_df


def plot_openpmd(args):
    z_slice_value = float(args["--z-slice"])

    verbose = bool(args["--verbose"])
    openpmd_format = (args["--openpmd-format"])

    data_file = args["<data-file>"]
    fname = f"{data_file}.it%08T{openpmd_format}"

    thorn_name = args["<thorn-name>"]
    group_name = args["<group-name>"]
    var_name = args["<var-name>"]

    patch = args["<num-patches>"].zfill(2)
    int_patch = int(args["<num-patches>"])
    patch_range = range(int_patch + 1)

    level = args["--refinement-level"].zfill(2)
    iteration_index = int(args["<iteration>"])
    iteration_index_pad = args["<iteration>"].zfill(8)

    save_file = bool(args["--save"])
    diverging = bool(args["--diverging"])

    mesh_name = f"{thorn_name}_{group_name}_patch{patch}_lev{level}"
    gf_name = f"{thorn_name}_{var_name}"

    # Get dataframes for each patch in parallel
    with multiprocessing.Pool(processes=(int_patch + 1)) as pool:
        futures = [
            pool.apply_async(
                merge_multipatch_dataframes,
                [
                    verbose,
                    fname,
                    iteration_index,
                    str(p).zfill(2),
                    level,
                    mesh_name,
                    gf_name,
                    z_slice_value
                ]
            ) for p in patch_range
        ]

        plt.close("all")
        for future in futures:
            df = future.get()
            if diverging == True:
                plt.tricontourf(
                    df["coordinatesx_vcoordx"],
                    df["coordinatesx_vcoordy"],
                    df[gf_name],
                    cmap="seismic",
                    levels=np.linspace(-1, 1, 101)
                )
            else:
                plt.tricontourf(
                    df["coordinatesx_vcoordx"],
                    df["coordinatesx_vcoordy"],
                    df[gf_name],
                    levels=100,
                )

        cb = plt.colorbar()
        cb.ax.set_ylabel(var_name)

        plt.tight_layout()

        if not save_file:
            plt.show()
        else:
            plt.savefig(f"{var_name}_it{iteration_index_pad}.png")


def main(args):
    # MPL settings.
    mpl.rcParams['mathtext.fontset'] = 'cm'
    mpl.rcParams['font.family'] = 'Latin Modern Roman'

    if args["augment-tsv"]:
        augment_tsv(args["<data-file>"], args["<coord-file>"])
    elif args["plot-tsv"]:
        plot_tsv(args)
    elif args["plot-openpmd"]:
        plot_openpmd(args)


# Required in order to keep subprocesses from launching recursivelly
if __name__ == '__main__':
    arguments = docopt(__doc__, version="mpx 1.0")
    main(arguments)
