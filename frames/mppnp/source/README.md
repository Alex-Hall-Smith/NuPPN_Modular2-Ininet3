## Setting up the Makefile configuration files

`CODE/Make.local` is the Makefile for the local machine. See MAKE_LOCAL for
examples. It should contain the architecture of the local machine, the path to
the SE libraries, the MPI installation and, depending on what solver you will
be using, information about BLAS and LAPACK libraries.

The Makefiles in the `ARCH` directories contain compiler/platform/OS specific
compiler flags and settings. See the `README` in `ARCH` for the naming
convention. These should not need to be changed by the user.

### LAPACK and BLAS (KML and ACML dense solvers)

LAPACK_LIBS and LD_PATH_LAPACK are required for the Intel/ACML optimised dense
matrix inversion routines (solver option 3) and the Intel Direct Sparse Solver
(option 4). For this you should have Intel MKL or AMD's ACML available on your
machine. ACML is free but not open source, and one must buy a license for MKL.
For more info on how to set up ACML with NuGrid codes, see [this
page](http://nugridstars.org/work-packages/solver-package/requirements-for-acml-gfort).
Anyway, this is not the fastest solver option, see below.

### SUPERLU sparse solver

`SUPERLU_PATH` and `BLASLIB` are required for solver option 5 (superLU). This
is the fastest option we have and requires a local build of the superLU
libraries.  These are completely open source and are distributed with the
NuGrid solver package. One only needs a BLAS installation, which is linked at
compile time using the `BLASLIB` variable in the `Make.local` file.

##Compiling MPPNP

In the `CODE` directory of the mppnp frame (where this file lives), compiling
the code requires the file parameter.inc in order to know the sizes of certain
arrays. If this is not already here (which it should not be before you have
compiled), it is built from `parameter_frame.inc` (from this directory) and
parameter_physics.inc (from the physics `CODE` directory). These array sizes
are somewhat problem-specific and may need to be changed from run to run, for
example if the network size is changed. In the run directories these sizes are
set/detected in the `ppn_frame.input` and `ppn_physics.input` files.

* `make distclean`: clean out the frame/run directory as well as physics and
solver packages
* `make superclean`: clean out the superLU libraries
* `make clean`: only clean frame/run directory
* `make`: generate optimised production code
* `make debug`: generate code for debugging

## Running MPPNP:

See the `README` in the `RUN_TEMPLATE` directory on how to set up the run
directory and run the code.

##Testing h5 read/write:

It happens again and again that you compile and run (usually for the first time
on a new machine) and the h5 throws IO errors. You think you have not correctly
compiled the h5 (usually you have!). `se_read_test.f` is a little test program
that you can compile instead of `mppnp.f` and check a simple read and write.
You may have to edit that little test program a bit as it has, for example, a
hard wired absolute path name for the test se.h5 file.
