import argparse

__description = (
        "This estimates quantum tunneling properties of an particle"
        " passing through a 2D material sheet."
        "The model is described in"
        " J. Chem. Phys. 148 (2018) 204707.\n\n"
        "The necessary input is the mass of the tunneling particle"
        " and the total energy of the 2D material/particle system"
        " in function of the perpendicular distance from the 2D"
        " material to the particle."
)

__epilog = (
        "" #TODO
)

parser = argparse.ArgumentParser(
        prog = "qtp",
        description = __description,
        epilog = '\n'.join(__epilog),
        formatter_class = argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars = '@',
        )

parser.add_argument('datafile', metavar="<file>", type=str, 
                    help=("file containting the total energy (surface +"
                          " tunneling atom) in function of the `z` distance"
                          " between the atom and the 2D surface."
                          " The file must be in XY format, i.e., `z` values"
                          " in the first column and the corresponding"
                          " energies in the second column."
                          " Moreover, energies must be given in 'hartrees' and"
                          " distances in 'angstroms'."
                          " Note that `z` can be negative, depending on"
                          " whether the atom is above or below the surface."
                          " Furthermore, lines starting with `#` will be"
                          " ignored by the parser."))

parser.add_argument('-m', dest='pmass', metavar="<mass>", type=float, default=1,
                    help=("mass of the tunneling particle in 'daltons'"
                          " (1 Da = 1 g/mol). (Default: 1 Da)"))

parser.add_argument('-temp', dest='temp', metavar="<temperature>",
                    type=float, default=[300], nargs='*',
                    help=("temperature in 'kelvins'."
                          " It can be specified as a single value or as a"
                          " range defined by three numbers: lower limit,"
                          " upper limit, and step size."
                          " For example, `-temp 300 800 10` requests all"
                          " temperatures in the interval [300 K, 800 K)"
                          " at steps of 10 K."))

parser.add_argument('-zrange', dest='zrange_angst', metavar="<value>",
                    type=float, nargs=3, default=None, required=False,
                    help=("range of coordinates for which local"
                          " properties should be printed in the output."
                          " Like the temperature, the range of coordinates is"
                          " defined by 3 numbers: lower limit, upper limit,"
                          " and step size."
                          " E.g., `-z -0.5 2.0 0.2` requests all coordinates"
                          " from -0.5 bohr up to 2.0 bohr at intervals of"
                          " 0.2 bohr."
                          " If coordinates are not specified, the program"
                          " will use the values in the input data file."
                          " (N.B. Since the input data are converged into a"
                          " spline function, this can be used to produce a"
                          " smoother curve by specifying a finer z step to"
                          " interpolate between the input points.)"))

def main():
    return parser.parse_args()

if __name__ == "__main__":
    args = main()
    print(args)

