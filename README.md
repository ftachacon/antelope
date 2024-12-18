## Antelope: Attosecond electron dynamics in condensed-matter physics 

Antelope pursues the development of a general code to simulate the non-linear light emission responses from a solid while a strong and ultrashort femtosecond laser source interacts with the target. To this end, we solve the Semiconductor Bloch Equations (SBE) or the time-dependent density matrix (TDDM) of the system and compute the time-dependent inter- and intra-band currents of the medium.

You can download and clone typing via terminal:

$ git clone https://github.com/ftachacon/antelope.git



## Prerequisites

**System Requirements:**

antelope is known to work on GNU/Linux. However, it should work on any POSIX-compliant system.

**Dependencies:**

- **mpicxx** (ubuntu: sudo apt-get install mpi-default-dev)

- **python-tk** (ubuntu: sudo apt-get install python-tk)

- **matplotlib** (pip2 install matplotlib --user)

- **numpy** (pip2 install numpy --user)

- **[GNU Bash](https://www.gnu.org/software/bash/)**

- **[GNU Awk (gawk)](https://www.gnu.org/software/gawk/)** (version >= 4.0)

- **[Intel® Fortran Compiler](https://software.intel.com/en-us/fortran-compilers)** (version >= 14.0.3)<br>
  M3C has not been tested with any other compiler.
  
- **[Intel® Math Kernel Library (Intel® MKL)](https://software.intel.com/en-us/mkl)**<br>
  M3C has not been tested with any other math library.

- **[SciFT (Scientific Fortran Tools)](https://github.com/nfaguirrec/scift)**<br>

**Recommended Dependencies:**

These dependencies are optional, but strongly recommended for full functionality:

- **[libmsym](https://github.com/mcodev31/libmsym)**<br>
  libmsym is a C library dealing with point group symmetry in molecules.

- **[GAUSSIAN](http://gaussian.com/)**<br> (version >= g09)
  Gaussian is a general purpose computational chemistry software package.

- **[GAMESS](https://www.msg.chem.iastate.edu/gamess/)**<br>
  General Atomic and Molecular Electronic Structure System GAMESS(US) is a software for computational chemistry.

## Compiling M3C

Download the .zip file from this page and extract the files,
```
$ unzip M3C-master.zip 
Archive:  M3C-master.zip
9f1572142803f97705c5db16b2d018a9c853a658
   creating: M3C-master/
  inflating: M3C-master/LICENSE      
  inflating: M3C-master/LICENSE.jmol
...

$ mv M3C-master M3C
```
or clone the repository using git
```
$ git clone https://github.com/nfaguirrec/M3C.git
```
The following should be the content of the M3C directory if previous steps were successful:
```
$ cd M3C
$ ls
doc   doxyfile  LICENSE.jmol     LICENSE.scift  Makefile   src        utils
docs  LICENSE   LICENSE.libmsym  M3Cvars.sh     README.md  templates  VERSION
```

Enter in the M3C directory (`cd M3C`) and modify the Makefile file (`src/Makefile`). In particular, choose the right path to the scift library (`-I<PATH_TO_SCIFT>/src` and `-L<PATH_TO_SCIFT>/src`).

To build the code just type make inside the main directory as follows:
```
$ make
cd src; make; cd ..
make[1]: Entering directory '/scratch/nestor/M3C/src'
Building dependencies for AzizSlamanPotential.f90 ... OK
Building dependencies for Fragment.f90 ... OK
Building dependencies for FragmentsDB.f90 ... OK
...
Building MarkovChain.f90 (0:03.40)
Building NNLS.f90 (0:00.08)
Building M3CBR.f90 (0:01.53)
```

## Installing M3C

The basic environmental variables that M3C needs can be loaded just adding the following command anywhere in the ~/.bashrc file:

```
source <PATH_TO_M3C>/M3Cvars.sh
```

M3C is also able to obtain data from electronic structure calculations by interfacing with some standard quantum chemistry programs. To enable this option the following variables must be specified:

```
# GAMESS configuration
export M3C_GAMESS_HOME=<PATH_TO_GAMESS_INSTALLATION>
export M3C_GAMESS_SCRATCH=/scratch/$USER/gamess

# GAUSSIAN configuration
export M3C_GAUSSIAN_HOME=<PATH_TO_GAUSSIAN_INSTALLATION>
export M3C_GAUSSIAN_SCRATCH=/scratch/$USER/gaussian
```

## Authors

* Alexis Chacon ( ftachacon@gmail.com )
* Nestor F. Aguirre ( nfaguirrec@gmail.com )

## Citing


