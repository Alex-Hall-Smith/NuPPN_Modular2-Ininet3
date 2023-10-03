! Parameters, variables and switches for the frame
module frame_knobs
   use utils, only: r8, i4
   implicit none
   private

   integer(i4), public :: &
         iolevel, & ! degree of verbosity of the code
         iplot_flux_option = 0, &
         iabuini, &
         plist, &
         interpolation_type, & !< for trajectory. 1: linear; 3: akima
         istart = -1 !< ancient init flag

   real(r8), public :: &
      ye_initial, &      !< initial \f$ Y_\mathrm{e} \f$ for NSE calculations
      t_end              !< maximum simulation time (s); 1e99 by default

   character*80, public :: &
      ini_filename, & ! filename for initial abundances
      tppnp_input_filename !< file name for input file containing initial abundances and trajectories

   logical, public :: restart

   ! public routines
   public :: readframeinput

contains


   subroutine readframeinput(t9threshold, ythreshold, datdir)
         ! *** read ppn_frame.input
         use communication
         implicit none
         real(r8) :: t9threshold, ythreshold
         character*256 :: datdir

         namelist/ppn_frame/iolevel, t9threshold, ythreshold, datdir, iabuini, ye_initial, &
            tppnp_input_filename, interpolation_type, t_end

         ! *** defaults

         iolevel     = 3
         t9threshold = 0.008_r8
         ythreshold  = 1.d-99
         iabuini     = 101
         ye_initial  = 0.5d0
         interpolation_type = 3 !< akima
         t_end       = 1.e99_r8

         datdir      = './input_data'
         tppnp_input_filename = "tp.lsb"

         if (master) then
            open(334,file='ppn_frame.input')
            read(334,NML=ppn_frame)
            close(334)
         end if

         call broadcast(iolevel)
         call broadcast(t9threshold)
         call broadcast(ythreshold)
         call broadcast(iabuini)
         call broadcast(ye_initial)
         call broadcast(interpolation_type)
         call broadcast(t_end)

   end subroutine readframeinput

end module frame_knobs
