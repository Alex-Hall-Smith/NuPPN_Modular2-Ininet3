! parameters, switches and variables for the frame
module frame_knobs
   use utils
   integer :: iolevel ! how verbose should the terminal output be?
   integer :: iplot_flux_option ! 0: no print fluxes; 1: print fluxes in
   ! flux_*.DAT
   integer :: i_flux_integrated ! 0: no integrated fluxes; 1: integrated
   ! fluxes
   integer :: iabuini           ! choice of initial abundances
   integer :: plist             ! for limiting species written out
   character*80 :: ini_filename ! filename for initial abundances
   character*80 :: ini_filename2 ! filename for initial abundances
   character*80 :: cprefix      ! prefix of abundance output
   character*3 :: code_source
   real(r8), public :: ye_initial      !< initial \f$ Y_\mathrm{e} \f$ for NSE calculations

   ! ^_^ init flag
   integer :: istart = -1

end module frame_knobs
