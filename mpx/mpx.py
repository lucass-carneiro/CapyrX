"""mpx: CapyrX plot and validation utility.

Usage:
  mpx plot-openpmd-slice <data-file> <thorn-name> <group-name> <var-name> <iteration> <patches-data>
                         [--verbose] [--save] [--diverging] [--autorange] [--varmin=<min>] [--varmax=<max>]
                         [--slice-coord=<coord>] [--slice-val=<val>] [--refinement-level=<level>] [--openpmd-format=<format>]
                          
  mpx (-h | --help)
  mpx --version

Options:
  -h --help                   Show this screen.
  --version                   Show version.
  --verbose                   Log operations.
  --save                      Save plot to file.
  --diverging                 Use a diverging color map.
  --autorange                 Computes the colorbar range automatically. Requires --diverging. Ignores --varmin and --varmax.
  --varmin=<min>              Minimun value of the colorbar. Requires --diverging. [default: -1.0].
  --varmax=<max>              Maximum value of the colorbar. Requires --diverging. [default: 1.0].
  --slice-coord=<coord>       The coordinate direction to slice [default: z].
  --slice-val=<val>           The (global) value of the coordinate slice to plot [default: 0.0]
  --refinement-level=<level>  The refinement level to plot [default: 0].
  --openpmd-format=<format>   The underlying format of OpenPMD files [default: .bp5].
"""


import matplotlib as mpl
from docopt import docopt

import plot_openpmd as plt_opmd

import logging
logger = logging.getLogger(__name__)

import sys

def main(args):
    logging.basicConfig(
        format="[%(asctime)s] [PID: %(process)d] %(levelname)s: %(message)s",
        #filename="mpx.log",
        #filemode="w",
        stream=sys.stdout,
        level=logging.INFO,
    )

    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["font.family"] = "Latin Modern Roman"

    if args["plot-openpmd-slice"]:
        plt_opmd.plot_openpmd_slice(args)


# Required in order to keep subprocesses from launching recursivelly
if __name__ == "__main__":
    arguments = docopt(__doc__, version="mpx 1.0")
    main(arguments)
