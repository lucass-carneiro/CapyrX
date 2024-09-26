import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

import logging
logger = logging.getLogger(__name__)


def read_to_dataframes(coords_file: str, data_file: str):
    columns = [
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

    data_columns = [
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
        "u",
    ]

    df = pd.read_csv(coords_file, sep="\t", comment="#", names=columns)
    data_df = pd.read_csv(data_file, sep="\t", comment="#", names=data_columns)

    df["u"] = data_df["u"]

    return df


def filter_patch(df: pd.DataFrame, index: int):
    new_df = df[
        (df["patch"] == index) &

        ((df["x"] > -1.0) | np.isclose(df["x"], -1.0)) &
        ((df["x"] < 1.0) | np.isclose(df["x"], 1.0)) &

        ((df["y"] > -1.0) | np.isclose(df["y"], -1.0)) &
        ((df["y"] < 1.0) | np.isclose(df["y"], 1.0)) &

        ((df["z"] > -1.0) | np.isclose(df["z"], -1.0)) &
        ((df["z"] < 1.0) | np.isclose(df["z"], 1.0)) &

        np.isclose(df["vcoordz"], 0.0)
    ]
    return new_df


def plot_ascii(args):
    coords_file = args["<coords-file>"]
    data_file = args["<data-file>"]

    save_file = bool(args["--save"])

    diverging = bool(args["--diverging"])
    varmin = float(args["--varmin"])
    varmax = float(args["--varmax"])
    autorange = bool(args["--autorange"])

    data = read_to_dataframes(coords_file, data_file)

    data_list = [
        filter_patch(data, 0),
        # filter_patch(data, 1),
        # filter_patch(data, 2),
        # filter_patch(data, 3),
        # filter_patch(data, 4)
    ]

    if diverging == True:
        lvls = None

        if autorange:
            lvls = 100
        else:
            lvls = np.linspace(varmin, varmax, 101)
            for df in data_list:
                x = df["vcoordx"].to_numpy()
                y = df["vcoordy"].to_numpy()
                z = df["u"].to_numpy()

                triang = tri.Triangulation(x, y)
                # triang.set_mask(np.hypot(x[triang.triangles].mean(axis=1), y[triang.triangles].mean(axis=1)) < 1.0)

                plt.triplot(triang)
                plt.tricontourf(
                    triang,
                    z,
                    cmap="seismic",
                    levels=lvls
                )
    else:
        for df in data_list:
            plt.tricontourf(
                df["vcoordx"],
                df["vcoordy"],
                df["u"],
                levels=100,
            )

    plt.xlabel("x")
    plt.ylabel("y")

    # cb = plt.colorbar()
    # cb.ax.set_ylabel("u")

    plt.tight_layout()

    if not save_file:
        plt.show()
    else:
        plt.savefig(f"u.png")
