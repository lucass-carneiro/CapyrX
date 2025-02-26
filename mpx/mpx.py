"""mpx: CapyrX plot and validation utility.

Usage:
  mpx plot-openpmd-slice [options] <data-file> <thorn-name> <group-name> <var-name> <iteration> <patches-data>
  mpx conv-openpmd-slice [options] <coarse> <medium> <fine> <thorn-name> <group-name> <var-name> <iteration> <patches-data>
  mpx plot-ascii-slice [options] <coords-file> <data-file>
                          
  mpx (-h | --help)
  mpx --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --verbose                   Log operations.
  --save                      Save plot to file.
  --out-format=<fmt>          Format to use when saving images. [default: png]
  --out-dpi=<dpi>             Resolution of the output image, in dpi (dots per inch). Relevant only on raster formats, such as png. [default: 300]
  --diverging                 Use a diverging color map.
  --autorange                 Computes the colorbar range automatically. Requires --diverging. Ignores --varmin and --varmax.
  --varmin=<min>              Minimun value of the colorbar. Requires --diverging. [default: -1.0].
  --varmax=<max>              Maximum value of the colorbar. Requires --diverging. [default: 1.0].
  --slice-coord=<coord>       The coordinate direction to slice [default: z].
  --slice-val=<val>           The (global) value of the coordinate slice to plot [default: 0.0]
  --refinement-level=<level>  The refinement level to plot [default: 0].
  --openpmd-format=<format>   The underlying format of OpenPMD files [default: .bp5].
  --plot-tri                  Plots the underlying triangulation used by matplotlib.
  --mask-radius=<radius>      Masks a radius from the origin. Zero radius means no mask. Usefull in Thornburg06 [default: 0.0].
  --xmin=<min>                Min. x value in convergence plots.
  --xmax=<max>                Max. y value in convergence plots.
"""


import sys
import matplotlib as mpl
from docopt import docopt

import plot_openpmd as plt_opmd
import plot_ascii as plt_ascii
import conv

import logging
logger = logging.getLogger(__name__)


def main(args):
    logging.basicConfig(
        format="[%(asctime)s] [PID: %(process)d] %(levelname)s: %(message)s",
        # filename="mpx.log",
        # filemode="w",
        stream=sys.stdout,
        level=logging.INFO,
    )

    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["font.family"] = "Latin Modern Roman"

    if args["plot-openpmd-slice"]:
        plt_opmd.plot_openpmd_slice(args)
    elif args["plot-ascii-slice"]:
        plt_ascii.plot_ascii(args)
    elif args["conv-openpmd-slice"]:
        conv.conv_opmd_slice(args)

    # Required in order to keep subprocesses from launching recursivelly
if __name__ == "__main__":
    arguments = docopt(__doc__, version="mpx 1.0")
    main(arguments)
