# Running a new mppnp problem

1. Set up the Make.local file in the frames/mppnp/CODE directory for your
machine (the templates in the MAKE_LOCAL directory are a good starting point)

2. In the mppnp directoy, make a copy of the RUN_TEMPLATE

3. prepare your input files; the four input files that define the problem are:

    * ppn_frame.input:	everything related to the I/O, data handling etc
    * ppn_physics.input: determine nuclear physics input
    * ppn_solver.input: specify input parameters related to the numerical
solver
    * isotopedatabase.txt: isotopic data that the code needs to know about

4. In the new run directory compile the code (cleaning beforehand so that the
executable you produce is the correct one):

        make distclean
        make

5. To run, check first if there is a template job script or run script for your
machine in the jobscripts or scripts directories. If not, you can run the code
where you stand like so, where nprocs should be the number of processors you
want to use:

       $ mpirun -np nprocs ./mppnp.exe

6. Outputs and errors will be written according to the jobscript file if you
used one, or to the screen if you did not. Progress can be checked in the file
summaryinfo.dat.


#Plotting/Analysis: 

Use nugridse.py from [NuGridPy](https://github.com/NuGrid/NuGridPy)
