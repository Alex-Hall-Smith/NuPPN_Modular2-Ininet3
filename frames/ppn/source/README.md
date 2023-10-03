## Setting up the Makefile configuration files

In the Makefiles in `ARCH` directories we store only compiler/platform/OS
specific things. See the README in ARCH for the naming convention.

In `CODE/Make.local` we store things that are specific to a machine. 
There are several machine specific Make.local files (as for example
`Make.local.nubuntu-gfort-superLU-ACML-openBLAS`) in addition to a generic
`Make.local.sample`. These are stored in the `CODE/MAKE_LOCAL` directory. Copy whatever
suits your situation to `CODE/Make.local` and check that file before compiling.
Locations of things may have changed.



`USE_SUPERLU` and `BLASLIB` are required for solver option 5 (superLU). This is
the fastest option and requires a local build of the superLU libraries.
This happens automatically during the first compile if `USE_SUPERLU` is set to `YES`. `BLASLIB`
should point to the BLAS libraries. 


`LAPACK_LIBS` and `LD_PATH_LAPACK` are required for the solver option
(mat_solv_option in ppn_solver.input) 3 (and 4). For this you should have Intel
MKL or AMD's ACML available on your machine. ACML is free but not open source,
and one must buy a license for MKL. Anyway, this is slower than superLU.

Here are a few typical options:

```
# system blas (usually -L/usr/lib/libblas -lblas)   
# openblas libraries (-L/usr/lib -lopenblas)
#    these can be installed on debian/ubuntu by apt-get install libopenblas-dev
#    or you can check them out from utils/lib and point the linker to there
# intel MKL (-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack95_lp64 -lmkl_core
# -liomp5 -lpthread)
# ACML (-L/opt/acml-4-4-0/gfortran64/lib -lacml -lacml_mv)
# etc ...
```
