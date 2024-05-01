"""mpx: CapyrX plot and validation utility.

Usage:
  mpx augment-tsv <data-file> <coord-file>
  mpx plot-tsv <augmented-data> [--rhs] [--save] [--diverging]
  mpx plot-grid-tsv <coordinates-tsv-file> [--save]
  mpx plot-openpmd <data-file> <thorn-name> <group-name> <var-name> <iteration> <num-patches> [--refinement-level=<level>] [--z-slice=<zval>] [--openpmd-format=<format>] [--save] [--diverging] [--varmin=<min>] [--varmax=<max>] [--autorange] [--verbose]
  mpx (-h | --help)
  mpx --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --rhs                       Plot RHS variables.
  --save                      Save plot to file.
  --verbose                   Show OpenPMD file information.
  --diverging                 Use a diverging color map.
  --varmin=<min>              Minimun value of the colorbar. Requires --diverging. [default: -1.0].
  --varmax=<max>              Maximum value of the colorbar. Requires --diverging. [default: 1.0].
  --autorange                 Computes the colorbar range automatically. Requires --diverging. Ignores --varmin and --varmax.
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


def plot_grid_tsv(args):
    # Arguments
    save_file = arguments["--save"]
    coord_file = args["<coordinates-tsv-file>"]

    print("Plotting", coord_file)

    # Common file vars
    vars = [
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
        "vcoordx",
        "vcoordy",
        "vcoordz"
    ]

    print("Reading data")
    data = pd.read_csv(coord_file, sep=r"\s+",
                       names=vars, comment="#")

    print("Creating data plot")

    plt.close("all")

    # Plot patches with different colors and markers
    markers = [".", "x", "1", "s", "d", "v", "o"]

    for i in range(0, data["patch"].iloc[-1] + 1):
        patch_sliced_data = data.loc[
            (data["vcoordz"] - 0.0 < 1.0e-14) &
            (data["patch"] == i)
        ]

        plt.scatter(
            patch_sliced_data["vcoordx"],
            patch_sliced_data["vcoordy"],
            marker=markers[i],
            label="Patch {}".format(i)
        )

    plt.tight_layout()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.legend()

    if not save_file:
        plt.show()
    else:
        plt.savefig("grid.pdf")

    plt.close("all")


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
    save_file = arguments["--save"]
    diverging = arguments["--diverging"]

    # Datafile vars, for pandas
    if rhs:
        vars = vars1 + ["u_rhs", "rho_rhs"] + vars2
    else:
        vars = vars1 + ["u", "rho"] + vars2

    print("Reading data")
    data = pd.read_csv(augmented_state_data,
                       sep=r"\s+", names=vars, comment="#")

    # Z slice data
    print("Filtering data for the z = 0 slice")

    sliced_data = data.loc[
        np.isclose(data["vcoordz"], -0.1) &
        (np.abs(data["x"]) <= 1.0) &
        (np.abs(data["y"]) <= 1.0) &
        (np.abs(data["z"]) <= 1.0)
    ]

    # RHS error data
    if rhs:
        print("Creating error plot")

        expected = -3.0 * np.pi**2 * np.cos(np.pi * sliced_data["vcoordx"]) * np.cos(
            np.pi * sliced_data["vcoordy"])
        error = np.abs(expected - sliced_data["u"])

        max_error = np.amax(error)
        error = error / max_error

        make_plot(sliced_data["vcoordx"], sliced_data["vcoordy"],
                  error, "u" + "_err", save_file, diverging)

    print("Creating data plot")
    make_plot(sliced_data["vcoordx"], sliced_data["vcoordy"],
              sliced_data["u"], "u", save_file, diverging)


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

    all_data = []

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
        for k in chunck_k_range:
            for j in chunck_j_range:
                for i in chunck_i_range:
                    all_data.append([
                        i,
                        j,
                        k,
                        x0 + i * dx,
                        z0 + k * dz,
                        y0 + j * dy,
                        data[k, j, i]
                    ])

    if verbose:
        print("Merging dataframes")

    merged_df = pd.DataFrame(
        all_data,
        columns=["i", "j", "k", "x", "y", "z", gf_name]
    )

    merged_df = merged_df.sort_values(["i", "j", "k"]).reset_index(drop=True)
    
    return merged_df


def merge_multipatch_dataframes(verbose, fname, iteration_index, patch, level, thorn_name, group_name, gf_name, z_slice_value):
    # Coordiantes
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

    # Main data
    mesh_name = f"{thorn_name}_{group_name}_patch{patch}_lev{level}"

    if verbose:
        print("Loading main data", mesh_name, "/", gf_name)

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
    filtered_df = merged_df[np.isclose(merged_df["coordinatesx_vcoordz"], z_slice_value)]
    
    print(
        merged_df[
            np.isclose(merged_df["coordinatesx_vcoordz"], z_slice_value)
        ]
    )

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
    num_patches = int(args["<num-patches>"])

    #patches = ["{:02d}".format(i) for i in range(num_patches)]
    patches = ["01"]

    level = args["--refinement-level"].zfill(2)
    iteration_index = int(args["<iteration>"])
    iteration_index_pad = args["<iteration>"].zfill(8)

    save_file = bool(args["--save"])

    diverging = bool(args["--diverging"])
    varmin = float(args["--varmin"])
    varmax = float(args["--varmax"])
    autorange = bool(args["--autorange"])

    gf_name = f"{thorn_name}_{var_name}"

    # Get dataframes for each patch in parallel
    with multiprocessing.Pool(processes=len(patches)) as pool:
        results = [
            pool.apply_async(
                merge_multipatch_dataframes,
                [
                    verbose,
                    fname,
                    iteration_index,
                    p,
                    level,
                    thorn_name,
                    group_name,
                    gf_name,
                    z_slice_value
                ],
            )
            for p in patches
        ]

        plt.close("all")

        for result in results:
            df = result.get()

            if diverging == True:
                lvls = None

                if autorange:
                    lvls = 100
                else:
                    lvls = np.linspace(varmin, varmax, 101)

                plt.tricontourf(
                    df["coordinatesx_vcoordx"],
                    df["coordinatesx_vcoordy"],
                    df[gf_name],
                    cmap="seismic",
                    levels=lvls
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
    elif args["plot-grid-tsv"]:
        plot_grid_tsv(args)


# Required in order to keep subprocesses from launching recursivelly
if __name__ == '__main__':
    arguments = docopt(__doc__, version="mpx 1.0")
    main(arguments)
