# Setting up the Makefile configuration files:

1. copy an existing example file from `MAKE.LOCAL` to a file called
`Make.local` in the CODE directory
2. set the various variables for your system
3. make sure that the frame, physics and solver are clean by running `make
distclean`
4. compile by running `make`

LAPACK_LIBS and LD_PATH_LAPACK are required for the solver option
(mat_solv_option in ppn_solver.input) 3 (and 4). For this you should have Intel
MKL or AMD's ACML available on your machine. ACML is free but not open source,
and one must buy a license for MKL. Anyway, this is not the fastest solver
option, see below.

SUPERLU_PATH and BLASLIB are required for solver option 5 (superLU). This is
the fastest option we have and requires a local build of the superLU libraries.
These are completely open source and straightforward to setup. There are
separate instructions for this here: TODO: ADD INSTRUCTIONS; SOMEONE REMIND
SAM.

## other compilation options

`make distclean`: clean out the frame/run directory as well as physics and
solver packages  
`make superclean`: clean out the superLU libraries  
`make clean`: only clean frame/run directory  
`make`:	generate optimised production code  
`make debug`: generate non-optimised code for debugging  

# Running TPPNP:

cp RUN_TEMPLATE to some work directory name of your choice on the same
level, modify the *.input files for your problem, make the executable:

    $ cp -r RUN_TEMPLATE run1
    $ cd run1
    $ vim ppn_frame.input
    # ..... edit the input files .....
    $ make distclean
    $ make

To run, check first if there is a template job script or run script for your
machine in the jabscripts or scripts directories. If not, you can run the code
where you stand like so:

       $ mpirun -np <nprocs> ./tppnp

where <nprocs> should be the number of processors you want to use.

Outputs and errors will be written according to the jobscript file if you used
one, or to the screen if you did not.

