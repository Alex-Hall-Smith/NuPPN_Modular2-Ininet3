!> parameters, switches and variables for the frame
module frame_knobs
   use utils
   integer :: iolevel                !< how verbose should the terminal output be?
   integer :: iplot_flux_option      !< 0: no print fluxes; 1: print fluxes in flux_*
   integer :: i_flux_integrated      !< 0: no integrated fluxes; 1: integrated fluxes
   integer :: nsource                !< 0: hydrostatic burn; 1: external trajectory
   integer :: print_cycles           !< sparsity for cycle output
   integer :: terminal_cycles        !< sparsity for terminal output
   integer :: iabuini                !< choice of initial abundances
   integer :: interpolation_type     !< for trajectory. 1: linear; 3: akima
   character(len=80) :: ini_filename !< filename for initial abundances
   character(len=80) :: cprefix      !< prefix of abundance output
   real(r8), public :: ye_initial    !< initial \f$ Y_\mathrm{e} \f$ for NSE calculations
   logical :: print_rate             !< should ppn print rate vs temperature?
   integer :: which_rate             !< if printing rate, which one?
   ! neutrino properties
   real(r8), public :: &
      lnue, &        !< electron neutrino luminosity from PNS (erg/s)
      lnuebar, &     !< antielectron neutrino luminosity from PNS (erg/s)
      lnux, &        !< "x" (tau and mu) neutrino luminosity from PNS (erg/s)
      lnuxbar, &     !< anti-"x" neutrino luminosity from PNS (erg/s)
      enue, &        !< average energy of electron neutrino spectrum (MeV)
      enuebar, &     !< average energy of antielectron neutrino spectrum (MeV)
      enux, &        !< average energy of "x" neutrino spectrum (MeV)
      enuxbar, &     !< average energy of anti-"x" neutrino spectrum (MeV)
      nu_t, &        !< time (s) since neutrino photosphere shock breakout
      tau_nu, &      !< e-folding time for decay of neutrino flux from PNS
      r              !< radius of particle in explosion (cm)

   ! ^_^ init flag
   integer :: istart = -1

contains

   !> @brief read ppn_frame.input, the main input card
   subroutine readframeinput( t9, rho, dt, t_max, dt_max, dt_factor )
         use trajectories, only: trajectory_fn
         implicit none
         real(r8) :: t9, rho, dt, t_max, dt_max, dt_factor, t9_nw_ini, rho_nw_ini
         namelist/ppn_frame/ t9, rho, dt, t_max, dt_max, dt_factor, nsource, iabuini, ini_filename, iolevel, cprefix, &
               iplot_flux_option, i_flux_integrated, print_cycles, t9_nw_ini, rho_nw_ini, trajectory_fn, &
               ye_initial, interpolation_type, print_rate, which_rate, terminal_cycles, lnue, &
               lnuebar, lnux, lnuxbar, enue, enuebar, enux, enuxbar, nu_t, tau_nu

         t9                = 0.1_r8                  ! temperature
         rho               = 1.0e3_r8                ! density

         dt                = 1.0d-3                  ! initial timestep / yr
         t_max             = 1.0d+3                  ! max age of evolution / yr
         dt_max            = 1.0d+2                  ! max number dzeitj can be / yr
         dt_factor         = 1.5d0                   ! increase time step by this factor until dt_max is reached

         nsource           = 0                       ! source of outside trajectory: hydrostatic burn (0) or external trajectory(1)
         iabuini           = 1                       ! how to initialize
         ini_filename      = 'initial_abundance.dat'
         trajectory_fn     = 'trajectory.input'      ! trajectory file name
         iolevel           = 1                       ! controls input/output
         cprefix           = 'iso_massf'             ! prefix for abundance vector output (used
         iplot_flux_option = 0                       ! iplot_flux_option
         i_flux_integrated = 0                       ! if use flux for delta_t or integrated.

         print_cycles      = 20                      ! time cycles between printing abundance files
         terminal_cycles   = 20                      ! time cycles between printing terminal output
         t9_nw_ini         = 0.1d0                   ! temperature and density for which network
         rho_nw_ini        = 100.d0                  ! is initialized and na <sigma v> is reported in networksetup.txt

         ye_initial        = 0.5d0                   ! initial ye for NSE calculations
         interpolation_type = 0                      ! interpolation currently implemented but not
                                                     ! used (set to 0)

         print_rate        = .false.                 ! whether to write out reaction rate and stop the program
         which_rate        = -1                      ! which rate to write out, if print_rate

         lnue              = 1e52_r8                 ! (erg/s)
         lnuebar           = 1e52_r8                 ! (erg/s)
         lnux              = 1e52_r8                 ! (erg/s)
         lnuxbar           = 1e52_r8                 ! (erg/s)
         enue              = 3.724_r8*3.15_r8        ! (MeV) from Wilson and Mayle calculations
         enuebar           = 4.835_r8*3.25_r8        ! (MeV) from Wilson and Mayle calculations
         enux              = 8*3.15_r8               ! (MeV)
         enuxbar           = 8*3.25_r8               ! (MeV)
         nu_t              = 4.e-3_r8                ! (s) from Woosley et al. 1990
         tau_nu            = 3._r8                   ! (s)
         r                 = 4.793e6_r8              ! (cm)

         open(2,file='ppn_frame.input')
         read(2,NML=ppn_frame)
         close(2)
         return

   end subroutine readframeinput

end module frame_knobs
