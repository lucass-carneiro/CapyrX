import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import openpmd_utils as opmdu

import concurrent.futures

import logging
logger = logging.getLogger(__name__)


def conv_opmd_slice(args):
    slice_coord = str(args["--slice-coord"])
    slice_val = float(args["--slice-val"])

    verbose = bool(args["--verbose"])
    openpmd_format = (args["--openpmd-format"])

    file_c = args["<coarse>"]
    file_m = args["<medium>"]
    file_f = args["<fine>"]

    fname_c = f"{file_c}.it%08T{openpmd_format}"
    fname_m = f"{file_m}.it%08T{openpmd_format}"
    fname_f = f"{file_f}.it%08T{openpmd_format}"

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

    num_patches = len(patches)
    patch_dataframe_futures_c = [None] * num_patches
    patch_dataframe_futures_m = [None] * num_patches
    patch_dataframe_futures_f = [None] * num_patches

    # Get dataframes for each patch in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_patches) as executor:
        for i in range(0, num_patches):
            patch_dataframe_futures_c[i] = executor.submit(
                opmdu.merge_multipatch_dataframes,
                verbose,
                fname_c,
                iteration_index,
                patches[i],
                level,
                thorn_name,
                group_name,
                gf_name,
                slice_coord,
                slice_val
            )

            patch_dataframe_futures_m[i] = executor.submit(
                opmdu.merge_multipatch_dataframes,
                verbose,
                fname_m,
                iteration_index,
                patches[i],
                level,
                thorn_name,
                group_name,
                gf_name,
                slice_coord,
                slice_val
            )

            patch_dataframe_futures_f[i] = executor.submit(
                opmdu.merge_multipatch_dataframes,
                verbose,
                fname_f,
                iteration_index,
                patches[i],
                level,
                thorn_name,
                group_name,
                gf_name,
                slice_coord,
                slice_val
            )

    patch_dataframes_c = [f.result() for f in patch_dataframe_futures_c]
    patch_dataframes_m = [f.result() for f in patch_dataframe_futures_m]
    patch_dataframes_f = [f.result() for f in patch_dataframe_futures_f]

    opmdu.downsample_patch_dataframes(
        verbose, opmdu.Resolution.medium, slice_coord, patch_dataframes_m)

    opmdu.downsample_patch_dataframes(
        verbose, opmdu.Resolution.fine, slice_coord,  patch_dataframes_f)

    # autopep8: off
    if verbose:
        logger.info(f"Creating convergence plot for iteration {iteration_index}")
    # autopep8: on

    plt.close("all")

    for i in range(0, len(patch_dataframes_c)):
        df_c = patch_dataframes_c[i]
        df_m = patch_dataframes_m[i]
        df_f = patch_dataframes_f[i]

        x = df_c[global_x].to_numpy()
        y = df_c[global_y].to_numpy()

        gf_c = df_c[gf_name].to_numpy()
        gf_m = df_m[gf_name].to_numpy()
        gf_f = df_f[gf_name].to_numpy()

        conv = np.abs(np.abs(gf_m - gf_f) * 2.0**(4.0) - np.abs(gf_c - gf_m))

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
                lvls = 200
            else:
                lvls = np.linspace(varmin, varmax, 201)

            plt.tricontourf(
                triang,
                conv,
                cmap="seismic",
                levels=lvls
            )
        else:
            lvls = None

            if autorange:
                lvls = 200
            else:
                lvls = np.linspace(varmin, varmax, 201)

            plt.tricontourf(
                x,
                y,
                conv,
                levels=lvls
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

    if verbose:
        logger.info("Done")
