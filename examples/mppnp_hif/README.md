## about this example

The example calculates 2000 cycles of H ingestion
based on the AGB model RUN-49/strat-B from Herwig et al. 2011
which matches the abundance distribution of Sakurai's object (Fig. 9 in Herwig+11).

## specific insturctions for this example

* download the SE data and the restart file from the CADC with the 
`./setup.sh` script


# Additional details:

1. Create the se files by using write_se_cycles.py:

	This script first reads trajectory data
	from the e2D14.0077501.se.h5 file,
	then it creates a grid and interpolates the trajectory
	data onto the new grid. Finally, it writes out
	the data using the sewrite.py routine (on svn:pylib).

2. For a restart, write the restart file by using write_restart.py:

	This script first reads the initial abundances from the file 
        restart0077991.check, then it gets the se infos (e.g. header) 
        from the file e2D14.0077501.se.h5.
	After that, it creates the grid and writes the restart
	file using sewrite.py.
	
The hif on/off (1/0) switch is in the ppn_frame.input file.
HIF parameters can be changed in the additional file ppn_hif.input.


