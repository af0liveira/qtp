# Quantum Transport Properties (QTP)

Calculation of quantum transport properties of a particle traveling
perpendicularly through a 2D-material sheet. 

This program was written to study tunneling properties of H+ isotopes through
graphene monolayers, but should work with different materials as well. 
The implementation is based on the equations published in:  
Poltavsky, I.; Zheng, L.; Mortazavi, M.; Tkatchenko, A. 
"Quantum Tunneling of Thermal Protons through Pristine Graphene."
_The Journal of Chemical Physics_ 2018, 148, 204707
(DOI: 10.1063/1.5024317).


## Requirements

The list below reflects the resources used in the development of this
program. Although not tested, it is likely that other versions will work as
well.

* python 3.6.8
* numpy 1.16.2
* scipy 1.2.1


## Setting up the program

N.B. The instructions below assume that you are at the repository's base
directory.

An automatic installation procedure is not available for this program. 
Instead, the program needs to be downloaded and the scripts in the `bin`
directory should be changed to point to the correct python interpreter. 
To do so, edit the `./bin/qtp.sh` file as follows:

* `py3` must point to your Python 3 interpreter (might not require any changes
  if your default is `python3`);
* `scrpath` must point to the folder where you have the QTP modules -- i.e.,
  `./qtp`, which contains the `.py` files.

If you specify the variables above using full paths, you should be able to
copy the `./bin/qtp.sh` file to anywhere in your system and still be able
to run the program.


## Running the program

N.B. The instructions below assume that you are at the repository's base
directory.

As mentioned above, you can run the program by executing the `bin/qtp.sh`
file with the desirable command line arguments.
For example:
```
bin/qtp.sh potential.dat -m 3 -temp 100 300 10
```
will run the program for triton particles (m = 3 Da) and temperatures from
100 to 300 K at steps of 10 K, using the potential energy barrier specified
in the `potential.dat` file.

For more details, simply run:

`bin/qtp.sh -h`


## Input data

### Potential energy barrrier

The potential energy barrier is passed via a file containing two columns:

1. The distances _z_ from the maximum energy barrier
2. The corresponding potential energies _U_(_z_)

The value must be in angstroms and hartrees, respetively.

**Important:**
The input potential barrier will be converted into a spline function
following the steps:

1. The potential will be shifted so that the lowest energy will be zero.
2. The potential will be symmetrized around _z_=0; this is done by
   mirroring the points on the _U_(_z_) axis, thus ensuring
   _U_(-_z_) = _U_(_z_).

### Temperature ranges

Note that the program calculates activation energies and prefactors
numerically and, thus, at least 5 temperature values are required.
If the number of points is not enough, the program will skip this part of
the calculation.

---
