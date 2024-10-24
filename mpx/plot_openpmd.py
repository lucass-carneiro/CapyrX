import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import openpmd_utils as opmdu

import multiprocessing

# DEBUG
# import pandas as pd

import logging
logger = logging.getLogger(__name__)


def plot_openpmd_slice(args):
    slice_coord = str(args["--slice-coord"])
    slice_val = float(args["--slice-val"])

    verbose = bool(args["--verbose"])
    openpmd_format = (args["--openpmd-format"])

    data_file = args["<data-file>"]
    fname = f"{data_file}.it%08T{openpmd_format}"

    thorn_name = args["<thorn-name>"]
    group_name = args["<group-name>"]
    var_name = args["<var-name>"]
    patches = eval(args["<patches-data>"].encode("unicode_escape"))

    level = args["--refinement-level"].zfill(2)
    iteration_index = int(args["<iteration>"])
    iteration_index_pad = args["<iteration>"].zfill(8)

    save_file = bool(args["--save"])

    diverging = bool(args["--diverging"])
    varmin = float(args["--varmin"])
    varmax = float(args["--varmax"])
    autorange = bool(args["--autorange"])

    gf_name = f"{thorn_name}_{var_name}"

    plot_tri = bool(args["--plot-tri"])
    mask_radius = float(args["--mask-radius"])

    # Get the name of the global coordinate GFs used for the plots
    if slice_coord == "z":
        global_x = "coordinatesx_vcoordx"
        global_y = "coordinatesx_vcoordy"
    elif slice_coord == "y":
        global_x = "coordinatesx_vcoordx"
        global_y = "coordinatesx_vcoordz"
    elif slice_coord == "x":
        global_x = "coordinatesx_vcoordy"
        global_y = "coordinatesx_vcoordz"
    else:
        logger.error(f"Unrecognized slice coordinate {slice_coord}")
        raise RuntimeError(f"Unrecognized slice coordinate {slice_coord}")

    # Get dataframes for each patch in parallel
    with multiprocessing.Pool(processes=len(patches)) as pool:
        results = [
            pool.apply_async(
                opmdu.merge_multipatch_dataframes,
                [
                    verbose,
                    fname,
                    iteration_index,
                    p,
                    level,
                    thorn_name,
                    group_name,
                    gf_name,
                    slice_coord,
                    slice_val
                ],
            )
            for p in patches
        ]

        plt.close("all")

        for result in results:
            df = result.get()

            # DEBUG
            # pd.options.display.precision = 16
            # print(df[np.isclose(df["coordinatesx_vcoordz"], 1.0)])
            # exit()

            x = df[global_x].to_numpy()
            y = df[global_y].to_numpy()
            gf = df[gf_name].to_numpy()

            triang = tri.Triangulation(x, y)
            if (mask_radius > 0.0 and not np.isclose(mask_radius, 0.0)):
                triang.set_mask(
                    np.hypot(
                        x[triang.triangles].mean(axis=1),
                        y[triang.triangles].mean(axis=1)
                    ) < mask_radius
                )

            if diverging == True:
                lvls = None

                if autorange:
                    lvls = 100
                else:
                    lvls = np.linspace(varmin, varmax, 101)
                plt.tricontourf(
                    triang,
                    gf,
                    cmap="seismic",
                    levels=lvls
                )
            else:
                lvls = None

                if autorange:
                    lvls = 100
                else:
                    lvls = np.linspace(varmin, varmax, 101)

                plt.tricontourf(
                    df[global_x],
                    df[global_y],
                    df[gf_name],
                    levels=lvls,
                )

            if (plot_tri):
                plt.triplot(triang)

        plt.xlabel(global_x)
        plt.ylabel(global_y)

        cb = plt.colorbar()
        cb.ax.set_ylabel(var_name)

        plt.tight_layout()

        if not save_file:
            plt.show()
        else:
            plt.savefig(f"{var_name}_it{iteration_index_pad}.png")
