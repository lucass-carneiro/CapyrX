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

    varmin = float(args["--varmin"])
    varmax = float(args["--varmax"])
    xmin = args["--xmin"]
    xmax = args["--xmax"]
    autorange = bool(args["--autorange"])

    gf_name = f"{thorn_name}_{var_name}"

    # Get the name of the global coordinate GFs used for the plots
    if slice_coord == "yz":
        slice_coord = ["y", "z"]
        global_x = "coordinatesx_vcoordx"
    elif slice_coord == "xz":
        slice_coord = ["x", "z"]
        global_x = "coordinatesx_vcoordy"
    elif slice_coord == "xy":
        slice_coord = ["x", "y"]
        global_x = "coordinatesx_vcoordz"
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

        gf_c = df_c[gf_name].to_numpy()
        gf_m = df_m[gf_name].to_numpy()
        gf_f = df_f[gf_name].to_numpy()

        np.seterr(divide='ignore', invalid='ignore')
        conv = np.abs((gf_c - gf_m) / (gf_m - gf_f))
        conv = np.log2(conv)

        plt.scatter(
            x,
            conv
        )

    plt.xlabel(global_x)
    plt.ylabel(var_name)

    if not autorange:
        plt.ylim(varmin, varmax)

    if xmin is not None and xmax is not None:
        plt.xlim(float(xmin), float(xmax))

    plt.tight_layout()

    if not save_file:
        plt.show()
    else:
        plt.savefig(f"conv_{var_name}_it{iteration_index_pad}.pdf")

    if verbose:
        logger.info("Done")
