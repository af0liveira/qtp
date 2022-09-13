import os
import sys

import datetime

import numpy as np

import warnings

from scipy import constants

sys.path.insert(0, os.path.dirname(__file__))
from cmdline import parser
from peb import PEB
from traco import TransCoeff
from pflux import PFlux
from arrhenius import Arrhenius
sys.path.pop(0)

dalton2me = constants.atomic_mass/constants.m_e
kelvin2au = 1 / (constants.value('Hartree energy')/constants.k)
angst2bohr = constants.angstrom/constants.value('Bohr radius')

def main():

    title = "Quantum Transport Properties v.3 (QTP3)"
    print(len(title) * "=")
    print(title)
    print(len(title) * "=")
    print("Timestamp:", datetime.datetime.now(), end=2*"\n", flush=True)
    print("N.B. Unspecified units imply atomic units (hartree, bohr, etc.)", 
          end=2*"\n", flush=True)


    #-------------------- 
    # Parse input data
    #-------------------- 
    # N.B., suffixes are used to indicate values _not_ in atomic units.

    args = parser.parse_args()

    pmass_Da = args.pmass  # particle's mass in daltons
    pmass = pmass_Da * dalton2me    # particles mass in electron mass

    print(f"Particle's mass (Da):  {pmass_Da}",
          f"Particle's mass (m_e): {pmass}",
          sep="\n", end=2*"\n", flush=True)

    # temperatures in kelvin
    if len(args.temp) == 1:
        temps_K = args.temp
    elif len(args.temp) == 3:
        temps_K = list(np.arange(*args.temp))
    else:
        raise AssertionError(f"Invalid temperature input:"
                             f" {' '.join(args.temp)}")

    # z and U(z) values
    with open(args.datafile, 'r') as fp:
        values = [float(v) for s in fp.readlines()
                    if len(s.strip()) > 0 and s.strip()[0] != '#'
                    for v in s.strip().split()]

    zvals_angst = values[0::2]   # values of z in angstrom
    zvals = [z*angst2bohr for z in zvals_angst]

    Uvals = values[1::2]   # values of U(z) in hartree

    # instantiate PEB, TransCoeff, and PFlux
    peb = PEB(zvals, Uvals)
    traco = TransCoeff(peb, pmass)
    pflux = PFlux(traco)

    #-------------------- 
    # Print U(z) data
    #-------------------- 

    # Define coordinates for which data should be printed in the output
    if args.zrange_angst is not None:
        zcoords = list(np.arange(*args.zrange_angst)*angst2bohr)
        print("Info for z/angstrom from {x[0]:} to {x[1]:} (step={x[2]:})"\
                  .format(x=args.zrange_angst), flush=True)
    else:
        zcoords = zvals
        print(f"Info for all points in `{args.datafile}`", flush=True)

    head = "  ".join(
        [
            '{:^12s}'.format("z/angstrom"),
            '{:^12s}'.format("z/bohr"),
            '{:^12s}'.format("U(z)/hartree"),
            '{:^12s}'.format("dU/dz"),
            '{:^12s}'.format("T(U(z))"),
            '{:^12s}'.format("ln T(U(z))"),
        ]
    )
    print(len(head)*"-", flush=True)
    print(head, flush=True)
    print(len(head)*"-", flush=True)

    # zlist = [zvals[0] - (zvals[1]-zvals[0])]
    # zlist += zvals
    # zlist += [zvals[-1] + (zvals[-1]-zvals[-2])]
    # for z in zlist:
    for z in zcoords:
        ln_T, T = traco(z)
        row = "  ".join(
            [
                '{:>12.6f}'.format(z/angst2bohr),
                '{:>12.6f}'.format(z),
                '{:>12.6f}'.format(peb(z)),
                '{:>12.6f}'.format(peb(z, der=1)),
                '{:>12.6e}'.format(T),
                '{:>12.6f}'.format(ln_T),
            ]
        )
        print(row, flush=True)

    print(len(head)*"-", flush=True)
    print("Timestamp:", datetime.datetime.now(), end=2*"\n", flush=True)

    #----------------------------
    # Calculate properties of 1/T
    #----------------------------

    head = "  ".join(
        [
            '{:^10s}'.format("Temp/K"),
            '{:^10s}'.format("beta/(1/K)"),
            '{:^10s}'.format("beta/a.u."),
            '{:^14s}'.format("k (classic)"),
            '{:^14s}'.format("k (tunnel)"),
            '{:^14s}'.format("k (total)"),
            '{:^14s}'.format("flux (classic)"),
            '{:^14s}'.format("flux (tunnel)"),
            '{:^14s}'.format("flux (total)"),
        ]
    )
    print(len(head)*"-", flush=True)
    print(head, flush=True)
    print(len(head)*"-", flush=True)

    # these four lists will be used by the arrhenius module
    betas = []
    ln_kcs = []
    ln_kqs = []
    ln_ktots = []

    for t in temps_K:
        beta = 1/(t*kelvin2au)
        j_c, j_q, j_tot = pflux(beta)
        j2k = np.sqrt(2*np.pi*pmass*beta)
        k_c = j_c * j2k
        k_q = j_q * j2k
        k_tot = j_tot * j2k
        row = "  ".join(
            [
                '{:>10.2f}'.format(t),
                '{:>10.2e}'.format(1/t),
                '{:>10.2f}'.format(beta),
                '{:>14.6e}'.format(k_c),
                '{:>14.6e}'.format(k_q),
                '{:>14.6e}'.format(k_tot),
                '{:>14.6e}'.format(j_c),
                '{:>14.6e}'.format(j_q),
                '{:>14.6e}'.format(j_tot),
            ]
        )
        print(row, flush=True)
        betas.append(beta)
        ln_kcs.append(np.log(k_c))
        ln_kqs.append(np.log(k_q))
        ln_ktots.append(np.log(k_tot))
    print(len(head)*"-", flush=True)
    print("Timestamp:", datetime.datetime.now(), end=2*"\n", flush=True)

    #------------------------------
    # Calculate activation energies
    #------------------------------

    if len(temps_K) < 4:
        warnings.warn("Skipping activation energies."
                      "Less than 04 data points available.")
        sys.exit()

    arrhenius_c = Arrhenius(betas, ln_kcs)
    arrhenius_q = Arrhenius(betas, ln_kqs) 
    arrhenius_tot = Arrhenius(betas, ln_ktots)

    head = "  ".join(
        [
            '{:^10s}'.format('Temp/K'),
            '{:^10s}'.format('beta/(1/K)'),
            '{:^10s}'.format('beta/a.u.'),
            '{:^15s}'.format('Eact (classic)'),
            '{:^15s}'.format('Eact (tunnel)'),
            '{:^15s}'.format('Eact (total)'),
            '{:^15s}'.format('Coeff (classic)'),
            '{:^15s}'.format('Coeff (tunnel)'),
            '{:^15s}'.format('Coeff (total)'),
        ]
    )
    print(len(head)*"-", flush=True)
    print(head, flush=True)
    print(len(head)*"-", flush=True)
        
    for t, beta in zip(temps_K, betas):
        eact_c, coeff_c = arrhenius_c(beta)
        eact_q, coeff_q = arrhenius_q(beta)
        eact_tot, coeff_tot = arrhenius_tot(beta)
        row = "  ".join(
            [
                '{:>10.2f}'.format(t),
                '{:>10.2e}'.format(1/t),
                '{:>10.2f}'.format(beta),
                '{:>15.8f}'.format(eact_c),
                '{:>15.8f}'.format(eact_q),
                '{:>15.8f}'.format(eact_tot),
                '{:>15.6e}'.format(coeff_c),
                '{:>15.6e}'.format(coeff_q),
                '{:>15.6e}'.format(coeff_tot),
            ]
        )
        print(row, flush=True)
    print(len(head)*"-", flush=True)
    print("Timestamp:", datetime.datetime.now(), end=2*"\n", flush=True)


if __name__ == '__main__':
    main()
