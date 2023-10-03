#Prepare run by downloading input files from CADC

#get restart file
wget --content-disposition -q http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/ppd_wind/M5.00Z2.0e-02.standard/H5_restart/M5.00Z0.020.0001001.restart.h5?view=data
mkdir H5_restart
mv -v M5.00Z0.020.0001001.restart.h5 H5_restart/

#get SE input
wget -q --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M5.00Z2.0e-02/M5.00Z2.0e-02/M5.00Z0.020.idx?view=data
mkdir SE
mv -v M5.00Z0.020.idx SE/
wget -q --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M5.00Z2.0e-02/M5.00Z2.0e-02/M5.00Z0.020.0001001.se.h5?view=data
mv -v M5.00Z0.020.0001001.se.h5 SE/
wget -q --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/set1/set1.2/see_wind/M5.00Z2.0e-02/M5.00Z2.0e-02/M5.00Z0.020.0002001.se.h5?view=data
mv -v M5.00Z0.020.0002001.se.h5 SE/

