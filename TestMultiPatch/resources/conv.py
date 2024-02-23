import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_file(cols, file):
    df = pd.read_csv(file + ".txt", sep=" ", names=cols, header=None)
    df["error"] = np.abs(df["actual"] - df["expected"])
    filtered_df = df.loc[(df["patch"] == 0) & (df["j"] == 0) & (df["k"] == 0)]

    plt.close("all")
    plt.scatter(filtered_df["x"], filtered_df["error"])
    plt.tight_layout()
    plt.savefig(file + ".pdf")


def plot_conv(cols, files):
    dfs = [pd.read_csv(file + ".txt", sep=" ", names=cols,
                       header=None) for file in files]

    for df in dfs:
        df["error"] = np.abs(df["actual"] - df["expected"])

    filtered_dfs = [df.loc[(df["patch"] == 0) & (df["j"] == 0)
                           & (df["k"] == 0)] for df in dfs]
    
    coarse_df = filtered_dfs[0].iloc[[1,2,3,4,5,6]]
    medium_df = filtered_dfs[1].iloc[[1,3,5,7,9,11]]
    fine_df = filtered_dfs[2].iloc[[1,5,9,13,17,21]]

    x = coarse_df["x"]
    c_m_m = coarse_df["actual"].to_numpy() - medium_df["actual"].to_numpy()
    m_m_f = medium_df["actual"].to_numpy() - fine_df["actual"].to_numpy()
    
    conv = np.log2(c_m_m / m_m_f)

    plt.close("all")
    plt.scatter(x, coarse_df["error"], label="coarse")
    plt.scatter(x, medium_df["error"], label="medium")
    plt.scatter(x, fine_df["error"], label="fine")
    plt.legend()
    plt.tight_layout()
    plt.savefig("errors.pdf")

    plt.close("all")
    plt.scatter(x, conv)
    plt.tight_layout()
    plt.savefig("conv.pdf")
    


def main():
    cols = ["i", "j", "k", "patch", "x", "y", "z", "expected", "actual"]
    files = ["5", "10", "20"]

    for file in files:
        plot_file(cols, file)

    plot_conv(cols, files)


main()
