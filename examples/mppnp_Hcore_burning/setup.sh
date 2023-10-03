#To prepare for the example download the data from CADC

##prepare restart file
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/ppd_wind/M3.00Z2.0e-02.standard/H5_restart/M3.00Z0.020.0000001.restart.h5?view=data
#mkdir -v H5_restart
mv -v M3.00Z0.020.0000001.restart.h5 H5_restart/

##get SE input
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M3.00Z2.0e-02/M3.00Z2.0e-02/M3.00Z0.020.idx?view=data
mkdir -v SE
mv -v M3.00Z0.020.idx SE/
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M3.00Z2.0e-02/M3.00Z2.0e-02/M3.00Z0.020.0000001.se.h5?view=data
mv -v M3.00Z0.020.0000001.se.h5 SE/
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M3.00Z2.0e-02/M3.00Z2.0e-02/M3.00Z0.020.0001001.se.h5?view=data
mv -v M3.00Z0.020.0001001.se.h5 SE/

