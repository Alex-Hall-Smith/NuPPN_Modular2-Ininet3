      program se_read_test

      implicit none
      include 'parameter.inc'
      include 'mpif.h'
      
      include 'FSE.fi'

      integer FID1,FID2,mi,icountmod
      character*256 filein
      character*80 codev
      character*11 modname

      double precision mini,zini,rotini,overini,age_unit,mass_unit
      double precision radius_unit,rho_unit,temperature_unit
      double precision pressure_unit,velocity_unit,dcoeff_unit
      double precision dq(msl)

      icountmod = 1
      filein='/home/nugrid/data/M2.00Z0.010/M2.00Z0.010.0000001.se.h5'
      WRITE(6,*)"Reading from ", filein

      call FSE_OPEN(filein,FID1)

C *** reading in global parameters:
      call FSE_READ_SATTR(FID1,-1,"codev",codev)
      write(*,*) "Codev is = ",codev

      call FSE_READ_SATTR(FID1,-1,"modname",modname)
      call FSE_READ_DATTR(FID1,-1,"mini",mini)
      call FSE_READ_DATTR(FID1,-1,"zini",zini)
      call FSE_READ_DATTR(FID1,-1,"rotini",rotini)
      call FSE_READ_DATTR(FID1,-1,"overini",overini)

c Units for variables (with respect to CGS)
      call FSE_READ_DATTR(FID1, -1,"age_unit",age_unit)
      call FSE_READ_DATTR(FID1, -1, "mass_unit",mass_unit)
      call FSE_READ_DATTR(FID1, -1, "radius_unit",radius_unit)
      call FSE_READ_DATTR(FID1, -1, "rho_unit",rho_unit)
      call FSE_READ_DATTR(FID1, -1, "temperature_unit",temperature_unit)
c      call FSE_READ_DATTR(FID1, -1, "pressure_unit",pressure_unit)
c      call FSE_READ_DATTR(FID1, -1, "velocity_unit",velocity_unit)
      call FSE_READ_DATTR(FID1, -1, "dcoeff_unit",dcoeff_unit)      
      
      print *,"Dcoeff_unit = ", dcoeff_unit
      call FSE_READ_IATTR(FID1, icountmod, "shellnb", mi)
      print *,"number of shells = ", mi

      call FSE_READ_D(FID1, icountmod, mi, "delta_mass", dq)

      print *,"dq(1) = ", dq(1)

      call FSE_OPEN('test_out.h5',FID2)
      call FSE_WRITE_IATTR(FID2,-1,"numcodev",5)
      call FSE_WRITE_SATTR(FID2,-1,"codev","this code")
      call FSE_WRITE_DATTR(FID2,-1,"mini",2.987D0)
      call FSE_CLOSE(FID2)

      end
