import numpy as np
import pandas as pd
import openpmd_api as io

import concurrent.futures

import logging
logger = logging.getLogger(__name__)


def openpmd_to_dataframe(fname, verbose, iteration_index, mesh_name, gf_name):
    # Based on
    # https://gist.github.com/stevenrbrandt/ef3718d52680d8a4e6accbad198a0b42
    # with a few extra spices
    series = io.Series(fname, io.Access.read_only)

    if verbose:
        logger.info("Reading dataset")

    if iteration_index not in series.iterations:
        logger.error(f"Iteration {iteration_index} not available")
        raise RuntimeError(f"Iteration {iteration_index} not available")

    if iteration_index not in series.iterations:
        logger.error(f"Iteration {iteration_index} not available")
        raise RuntimeError(f"Iteration {iteration_index} not available")

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
        logger.info(
            # fmt: off
            f"Merging dataframes:\n  File: {fname}\n  Iteration: {iteration_index}\n  Mesh: {mesh_name}\n  Grid function: {gf_name}"
            # fmt: on
        )

    merged_df = pd.DataFrame(
        all_data,
        columns=["i", "j", "k", "x", "y", "z", gf_name]
    )

    merged_df = merged_df.sort_values(["i", "j", "k"]).reset_index(drop=True)
    return merged_df


def filter_multipatch_dataframes(verbose, patch, merged_df, slice_coord, slice_val):
    if verbose:
        logger.info(f"Filtering main data for {slice_coord} = {slice_val}")

    if (patch[1] == False):
        if slice_coord == "z":
            filtered_df = merged_df[np.isclose(merged_df["x"], slice_val)]
        elif slice_coord == "y":
            filtered_df = merged_df[np.isclose(merged_df["y"], slice_val)]
        elif slice_coord == "x":
            filtered_df = merged_df[np.isclose(merged_df["x"], slice_val)]
        else:
            logger.error(f"Unrecognized slice coordinate {slice_coord}")
            raise RuntimeError(f"Unrecognized slice coordinate {slice_coord}")
    else:
        if slice_coord == "z":
            filtered_df = merged_df[np.isclose(
                merged_df["coordinatesx_vcoordz"], slice_val)]
        elif slice_coord == "y":
            filtered_df = merged_df[np.isclose(
                merged_df["coordinatesx_vcoordy"], slice_val)]
        elif slice_coord == "x":
            filtered_df = merged_df[np.isclose(
                merged_df["coordinatesx_vcoordx"], slice_val)]
        else:
            logger.error(f"Unrecognized slice coordinate {slice_coord}")
            raise RuntimeError(f"Unrecognized slice coordinate {slice_coord}")

    return filtered_df


def merge_multipatch_dataframes(verbose, fname, iteration_index, patch, level, thorn_name, group_name, gf_name, slice_coord, slice_val):
    # Coordiantes
    if verbose:
        # autopep8: off
        logger.info(f"Loading data for iteration {iteration_index} patch {patch}")
        # autopep8: on

    mesh_name = f"{thorn_name}_{group_name}_patch{patch[0]}_lev{level}"

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        vcoordx_future = executor.submit(
            openpmd_to_dataframe,
            fname,
            verbose,
            iteration_index,
            f"coordinatesx_vertex_coords_patch{patch[0]}_lev{level}",
            "coordinatesx_vcoordx"
        )

        vcoordy_future = executor.submit(
            openpmd_to_dataframe,
            fname,
            verbose,
            iteration_index,
            f"coordinatesx_vertex_coords_patch{patch[0]}_lev{level}",
            "coordinatesx_vcoordy"
        )

        vcoordz_future = executor.submit(
            openpmd_to_dataframe,
            fname,
            verbose,
            iteration_index,
            f"coordinatesx_vertex_coords_patch{patch[0]}_lev{level}",
            "coordinatesx_vcoordz"
        )

        main_gf_future = executor.submit(
            openpmd_to_dataframe,
            fname,
            verbose,
            iteration_index,
            mesh_name,
            gf_name
        )

    vcoordx_df = vcoordx_future.result()
    vcoordy_df = vcoordy_future.result()
    vcoordz_df = vcoordz_future.result()
    main_gf_df = main_gf_future.result()

    if verbose:
        # autopep8: off
        logger.info(f"Merging data for iteration {iteration_index} patch {patch} ")
        # autopep8: on

    merge_cols = ["i", "j", "k", "x", "y", "z"]

    coords_df = pd.merge(vcoordx_df, vcoordy_df, on=merge_cols)
    coords_df = pd.merge(coords_df, vcoordz_df, on=merge_cols)
    merged_df = pd.merge(coords_df, main_gf_df, on=merge_cols)

    if verbose:
        # autopep8: off
        logger.info(f"Filtering ghost zones for iteration {iteration_index} patch {patch} ")
        # autopep8: on

    merged_df = merged_df[
        ((merged_df["x"] > -1.0) | np.isclose(merged_df["x"], -1.0)) &
        ((merged_df["x"] < 1.0) | np.isclose(merged_df["x"], 1.0)) &

        ((merged_df["y"] > -1.0) | np.isclose(merged_df["y"], -1.0)) &
        ((merged_df["y"] < 1.0) | np.isclose(merged_df["y"], 1.0)) &

        ((merged_df["z"] > -1.0) | np.isclose(merged_df["z"], -1.0)) &
        ((merged_df["z"] < 1.0) | np.isclose(merged_df["z"], 1.0))
    ]

    return filter_multipatch_dataframes(verbose, patch, merged_df, slice_coord, slice_val)
