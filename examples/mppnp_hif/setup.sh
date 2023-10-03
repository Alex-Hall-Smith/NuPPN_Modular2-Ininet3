#To download input data from the CADC for the run

#download sewrite.py and other packages from NuSE repo
wget https://raw.githubusercontent.com/NuGrid/NuSE/master/write_se/sewrite.py
wget https://raw.githubusercontent.com/NuGrid/NuSE/master/write_se/se.py
wget https://raw.githubusercontent.com/NuGrid/NuSE/master/write_se/h5T.py
#ascii_table is not found in regual nugridpy package so download here for now.
wget https://raw.githubusercontent.com/NuGrid/NuGridPy/master/ascii_table.py

#create restart file
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/projects/mppnp/examples/mppnp_hif/restart0077991.check?view=data
wget --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/projects/mppnp/examples/mppnp_hif/e2D14.0077501.se.h5?view=data
python  write_restart.py 

mkdir H5_restart
mv e2D14_hif* H5_restart/

#write SE input

mkdir SE
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/projects/mppnp/examples/mppnp_hif/e2D14_hif.idx?view=data
mv e2D14_hif.idx SE/
python write_se_cycles.py
mv e2D14_hif* SE/

#clean up
rm e2D14.0077501.se.h5
rm restart0077991.check
rm h5Preproc.txt
rm sewrite.py*
rm se.py*
rm h5T.py*

