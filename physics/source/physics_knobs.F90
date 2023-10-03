module physics_knobs
   use utils, only: r8
   use array_sizes, only: nsp
   use constants
   implicit none
   private

   integer, public :: ininet

   ! i_nse = 0 if t9 > nse_min_t9, the standard network is used without nse assumptions
   ! i_nse = 1 if t9 > nse_min_t9, nse equations are solved to give the abundances
   integer, public :: i_nse !< 0: don't use NSE, 1: use nse for t9 > nse_min_t9

   integer, public :: nse_option !< solver option for NSE: nse_torch (0) or nse_swj (1)
   integer, parameter, public :: NSE_TORCH = 0, NSE_SWJ = 1

   integer, public :: screen_option !< which screening to use: none(0), graboske(1) or chugunov(2)
   integer, parameter, public :: NO_SCREENING = 0, GRABOSKE = 1, CHUGUNOV = 2

   integer, public :: kadonis_interp ! 1 is linear, 3 is akima

   integer, public :: &
         nvcp, &
         nrcp, &
         nnn

   real(r8), public :: tbetamin
   real(r8), public :: decay_time ! time over which to decay composition if decay == .true.

   ! ^_^ switches
   ! whether or not to calculate the temperature and density derivatives of the reaction rates
   ! (will eventually be decided by the solver type rather than a hard switch).
   logical, parameter, public :: getderivs = .false.
   logical, public :: use_cache ! for reaclib
   logical, public :: detailed_balance ! get reverse rates via detailed balance
   logical, public :: weak_rates ! include weak rates in network?
   logical, public :: strong_rates ! include strong rates in network?
   logical, public :: decay ! is this a decay-only calculation?
   logical, public :: use_other_nuc ! whether to use reactions from the "other_nuc" module
   logical, public :: do_neutrinos ! whether to include neutrino reactions


   ! ^_^ index_reaclib: which reaclib to use
   !     0 = "Basel" Reaclib
   !     1 = JINA Reaclib
   !     2 = "Pure" Reaclib
   integer, public :: index_reaclib

   ! ^_^ CONTROLS:
   !     jbj mode is set in ppn_physics,input (see physics_knobs.F90 in the physics package)
   !     since JBJ16 rates cover the entire region of the isotopic chart that the Oda+ '94 rates cover, you can either:
   !        jbj_mode = 1 : don't use JBJ rates at all, just use Oda rates
   !        jbj_mode = 2 : use JBJ rates only where Oda rates are used, i.e. Oda rates are directly replaced
   !        jbj_mode = 3 : use full set of rates from Oda+ '94, and add remaining rates from JBJ where there is no overlap 
   !        jbj_mode = 4 : use JBJ rates everywhere provided, thus using none of the Oda rates (and extending beyond)
   integer, public :: jbj_mode

   !     NKK04 rate coverage overlaps with several other rate sources,
   !     so there are some options for how to merge the rates:
   !        nkk_mode = 1 : don't use NKK rates at all
   !        nkk_mode = 2 : use NKK rates only where LMP00 rates are available
   !        nkk_mode = 3 : use NKK rates when FFN, Oda94 and LMP are unavailable
   integer, public :: nkk_mode

   real(r8), public :: t9_nw_ini, nse_min_t9, ye_nw_ini, rho_nw_ini, yps_nw_ini(1,nsp)

   ! reaction rate multiplying factors
   integer(i4), parameter, public :: num_rate_factors = 10
   integer(i4), public :: rate_index(num_rate_factors)
   real(r8), public :: rate_factor(num_rate_factors)

   ! choice of treatment isomers
   ! 0 = do not include isomers
   ! 1 = include isomers
   integer, public :: isomer_choice 

   ! public routines
   public :: readphysicsinput

contains

   subroutine readphysicsinput()
         use communication

         NAMELIST/ppn_physics/ininet, i_nse, nvcp, nrcp, nnn, index_reaclib, tbetamin, jbj_mode, &
               nkk_mode, nse_min_t9, t9_nw_ini, ye_nw_ini, rho_nw_ini, use_cache, nse_option, &
               detailed_balance, screen_option, weak_rates, strong_rates, rate_index, rate_factor, &
               isomer_choice, decay, use_other_nuc, kadonis_interp, do_neutrinos, &
               decay_time
         integer nrcpold, fh
         common/old/nrcpold

         ininet              = 0
         i_nse               = 0
         nse_option          = 1
         nvcp                = 57
         nrcp                = 110
         nnn                 = 1107
         tbetamin            = 0.5_r8
         index_reaclib       = 2
         jbj_mode            = 1
         nkk_mode            = 1
         nse_min_t9          = 6._r8
         t9_nw_ini           = 7.8e-3_r8
         ye_nw_ini           = 0.5_r8
         rho_nw_ini          = 1._r8
         use_cache           = .false.
         yps_nw_ini(:,:)     = ZERO
         detailed_balance    = .true.
         screen_option       = CHUGUNOV
         weak_rates          = .true.
         strong_rates        = .true.
         rate_index(:)       = -1
         rate_factor(:)      = ONE
         use_other_nuc       = .false.
         isomer_choice       = 0
         decay               = .false.
         decay_time          = 1.e17 ! seconds
         kadonis_interp      = 3 ! 1 is linear, 3 is akima
         do_neutrinos        = .false.

         if ( master ) then
            open(2, file = "ppn_physics.input")
            read(2, nml = ppn_physics)
            ! if pure reaclib is used, then only JINAR can be used for now.
            ! Basel revision for now does not contain beta decays above Sn.
            if(ininet .eq. 4 .and. index_reaclib .eq. 0) then
               write(*,*) 'Basel reaclib does not include beta decays above Sn.'
               write(*,*) 'With ininet=4, use index_reaclib > 0.'
               stop
            end if
            if (detailed_balance .and. decay) then
               detailed_balance = .false.
               print *, "WARNING: detailed balance switched of for decays"
            end if
            ! I don't close file here because reading continues in rnetw2007
         end if

         nrcpold = nrcp

#ifndef PPN
         call broadcast(ininet)        ; call broadcast(i_nse)               ; call broadcast(nse_option)
         call broadcast(nvcp)          ; call broadcast(nrcp)                ; call broadcast(nnn)
         call broadcast(tbetamin)      ; call broadcast(index_reaclib)       ; call broadcast(jbj_mode)
         call broadcast(nkk_mode)      ; call broadcast(nse_min_t9)          ; call broadcast(t9_nw_ini)
         call broadcast(ye_nw_ini)     ; call broadcast(rho_nw_ini)          ; call broadcast(use_cache)
         call broadcast(yps_nw_ini)    ; call broadcast(detailed_balance)
         call broadcast(screen_option) ; call broadcast(weak_rates)          ; call broadcast(decay)
         call broadcast(rate_index)    ; call broadcast(use_other_nuc)       ; call broadcast(kadonis_interp)
         call broadcast(rate_factor)   ; call broadcast(do_neutrinos)
         call broadcast(isomer_choice) ; call broadcast(decay)               ; call broadcast(decay_time)
         call broadcast(strong_rates)
#endif

   end subroutine readphysicsinput


end module physics_knobs
