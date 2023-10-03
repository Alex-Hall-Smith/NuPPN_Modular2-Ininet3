#To download input data from the CADC for the run

#get restart file

wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.1/ppd_wind/M6.00Z1.0e-02/H5_restart/M6.00Z.0100.0004701.restart.h5?view=data

mv M6.00Z.0100.0004701.restart.h5 H5_restart/

#get SE input

mkdir SE
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.1/see_wind/M6.00Z1.0e-02/M6.00Z1.0e-02/M6.00Z.0100.idx?view=data
mv M6.00Z.0100.idx SE/

wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.1/see_wind/M6.00Z1.0e-02/M6.00Z1.0e-02/M6.00Z.0100.0004001.se.h5?view=data
mv M6.00Z.0100.0004001.se.h5 SE/

wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.1/see_wind/M6.00Z1.0e-02/M6.00Z1.0e-02/M6.00Z.0100.0005001.se.h5?view=data
mv M6.00Z.0100.0005001.se.h5 SE/


