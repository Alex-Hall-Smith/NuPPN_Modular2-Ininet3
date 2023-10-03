!
!  ooooooooo.   ooooooooo.   ooooo      ooo 
!  `888   `Y88. `888   `Y88. `888b.     `8` 
!   888   .d88`  888   .d88`  8 `88b.    8  
!   888ooo88P`   888ooo88P`   8   `88b.  8  
!   888          888          8     `88b.8  
!   888          888          8       `888  
!  o888o        o888o        o8o        `8  
!
!> @author Marco Pignatari, Falk Herwig, Samuel Jones and the NuGrid Collaboration
!
!> @brief
!> Single-zone nuclear reaction network code, based on ppn07-tdn9c.01.f, which was derived from nw6 by Falk Herwig.
!> Since developed by Marco Pignatari, Falk Herwig, Samuel Jones and the NuGrid Collaboration.
!

program ppn
   use array_sizes, only: nre, nsp
   use nuc_data, only: considerisotope, considerreaction, ispe, nuc_data_init, niso
   use rates, only: rates_init, v, v_rev, k1, k3, k5, k7
   use evaluate_rates, only: evaluate_all_rates
   use solver, only: solver_init, integrate_network
   use solver_diagnostics, only: ntime, nsubt, iter, titers, time_matrix_inv, &
         ntime_last
   use errors, only: MAX_SUBSTEPS, SMALL_DT
   use jac_rhs, only: makemap, nvar1
   use physics_knobs
   use frame_knobs
   use solver_knobs
   use utils, only: r8, i4
   use screening
   use reaclib, only: reaclib_init, reaclib_create_masks, reaclib_preprocessor
   use vital, only: vital_init, vital_rates_derivs
   use netgen, only: netgen_init
   use kadonis, only: kadonis_init
   use jbj16, only: jbj_init
   use fuller, only: fuller_init
   use nkk04, only: nkk_init
   use nse_torch, only: testnse
   use nse_swj, only: compute_nse, nse_init
   use alpha_decays
   use constants, only: avogadro, HALF, ONE, ZERO, yrs2sec
   use trajectories
   use output
   use abundances, only: set_initial_abundances
   use bader_deuflhard, only: bader_deuflhard_init
   use reaction_info, only: reaction_info_init, lab, labb, ilabb, rfac
   use physics, only: rnetw2007, rnetw2008
   use communication
   use reverse, only: reverse_init
   use other_nuc, only: other_nuc_init
   use networksetup, only: read_networksetup, nvnc1, nrnc1, write_networksetup
   use neutrinos, only: neutrinos_init, neutrino
   use decays, only: do_decay
   implicit none

   real(r8) :: yps(1,nsp), qi(nre), an(nsp), zn(nsp), t9, t9_0, t9_1, &
         rho, rho_0, rho_1, t_max_yr, dt, dt_yr, &
         dt_max_yr, agej, ye, oneyear, acsp(nsp), xndens, residual, dt_factor, tti1, tti2, &
         mu_p, mu_n, tstart, tsetup
   integer  :: nvar, nvrel, i, j, icountmod, ineut, iprint_true(nsp), ierr, iterations
   character(len=5) :: csp(nsp)

   character(len=5) :: zis_test(nsp), lab_replace(nre), labb_test(nre)
   real(r8) :: rfac_replace(nre)
   integer  :: ilabb_test(nre)
   logical  :: considerisotope_test(nsp), considerreaction_replace(nre)

   type(trajectory) :: t
   type(neutrino) :: nu

   common / cnetw / an, zn

   yps(1,:) = ZERO; icountmod = 0; ye = HALF
   tstart = wallclocktime()

   call print_splash
   print *, "ppn setting up network"

   call reaction_info_init()
   call solver_init()
   call vital_init()
   call fuller_init()
   call other_nuc_init()

   call readframeinput( t9, rho, dt_yr, t_max_yr, dt_max_yr, dt_factor )
   call readsolverinput()
   call readphysicsinput()

#if pIDX_RCLB == 3
   allocate(niso(0:i325dim,0:iCfdim,2))
#else
   allocate(niso(0:i282dim,0:iAtdim,2))
#endif

   ! only read tables one time
   istart = 0

   ! ^_^ initialisations - can these all go into rates_init or was there a gotcha?
   call screen_init()
   call rates_init()
   call reaclib_init()
   call netgen_init()
   call alpha_decays_init()
   call neutrinos_init(nu)
   call jbj_init()
   call nkk_init()
   call bader_deuflhard_init
   call kadonis_init()

   select case(ininet)
   case(3)
      ! network config read from networksetup.txt
      call read_networksetup(an, zn, t9_nw_ini, rho_nw_ini)
      !call rnetw2008( ye_nw_ini, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini, nu )

      ! *** initializing things to update and things to check for consistency between networksetup.txt and the physics package used...
      zis_test = zis
      considerisotope_test = considerisotope
      lab_replace = lab
      labb_test = labb
      ilabb_test = ilabb
      rfac_replace = rfac 
      considerreaction_replace = considerreaction
      call rnetw2007( ye_nw_ini, qi, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini ) 
      ! *** testing that nothing crazy was done by the editor...
      do i=1,nvnc1
       if (zis(i) /= zis_test(i) .or. considerisotope(i) .neqv. considerisotope_test(i)) then
        stop 'you modified the isotope names or removed them, you moron...'
       end if
      end do
      do i=1,nrnc1
       if (labb_test(i) /= labb(i) .or. ilabb(i) /= ilabb_test(i)) then
        print*,i,labb_test(i), labb(i), ilabb(i), ilabb_test(i)
        stop 'you are really stupid ... you modified labb or ilabb (@!~##!), or add by force reactions from the list?!'
       end if
      end do
      ! *** end of test. I am still alive, so I can replace now
      lab = lab_replace
      rfac = rfac_replace 
      considerreaction = considerreaction_replace
   case(4)
      ! pure reaclib
      stop 'rnetw_reaclib (ininet=4) in progress'
      ! call rnetw_reaclib( ye, qi, an, zn, nvar, nvrel, rho, t9, yps )
   case default
      ! set up network based on ppn_physics.input and isotopedatabase.txt
      call rnetw2007( ye_nw_ini, qi, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini )
   end select

   istart = 1

   call nuc_data_init()
   ! ^_^ make mapping from dense to sparse matrix format
   call makemap()
   if (detailed_balance) call reverse_init()

   ! ^_^ write networksetup file
   select case(ininet)
   case(3)
      ! write networksetup2
      ! stop "ppn frame: check souce code arguments to write_networksetup2"
      print *,'ininet=3 setting: all is ready to go'
   case(1)
      call write_networksetup(t9, rho, nvnc1, nrcp, nrnc1, an, zn)
      print *, 'stopping in rnetw2007:'
      print *, 'the file networksetup.txt has been created.'
      print *, 'you can make changes to the network by'
      print *, 'modifying this file. When you are happy with'
      print *, 'your changes, run the code again with ininet=3'
      print *, 'and your new networksetup.txt will be used.'
      stop
   case default
      call write_networksetup(t9, rho, nvnc1, nrcp, nrnc1, an, zn)
   end select

   ! ^_^ mask where reaclib rates are not included in the network, so that
   !     we do not waste time calculating them. This should really go into
   !     reaclib_init but with the current ordering in the way ppn
   !     sets up the network etc, this will have to be done later
   call reaclib_create_masks()
   if ( use_cache ) call reaclib_preprocessor

   ! iprint_true, gives the list of the .true. species
   ! here I am including and index to define what to print, only nvar species.
   iprint_true(:) = 0
   j = 0
   do i = 1, nsp
      if ( considerisotope(i) ) then
         j = j + 1
         iprint_true(j) = i
      end if
   end do

   ! initialize abundances
   call set_initial_abundances( yps )
   call nse_init()
   !> todo move this NSE initialisation into abundances module for ininet = 4
   if ( iabuini == 4 ) then
      print *, 'NSE_OPTION = ', nse_option
      select case( nse_option )
      case( NSE_TORCH )
         call testnse( t9, rho, ye_initial, nvar, yps, ierr )
      case( NSE_SWJ )
         call compute_nse( t9, rho, ye_initial, yps, mu_p, mu_n, iterations, ierr )
      case default
         stop 'NSE solver choice must be made'
      end select
   end if

   tti1 = wallclocktime()
   call check_mass_conservation(.true., yps, 1, residual)
   if (residual > 1.e-8_r8) call renorm_2d(yps, 1, residual, 'iniabund')
   call calculate_ye(yps(1,:), an, zn, ye, considerisotope)

   ! if decaying, do decay:
   if (decay) then
      call do_decay(nvar, yps, nvrel)
      call printonecycle(icountmod, dt_yr, agej, t9, rho, yps, acsp, csp, an, zn,.false.)
      stop "decays done"
   end if

   select case(nsource)
   case(0)
      print *, ' nsource = 0: hydrostatic burn (constant temperature and density)'
   case(1)
      print *, ' nsource = 1: run trajectory from file: ', trajectory_fn
      call ppn_trajectory_init(t, interpolation_type)
   case default
      stop 'nsource choice not supported; please choose 0 or 1'
   end select

   ! ^_^ initialization relic; TODO: can we get rid of this please.
   select case( nsource )
   case(0)
      agej = dt_yr
   case(1)
      agej = t%time(1) / yrs2sec
      t9   = t%temp(1)
      rho  = t%dens(1)
   end select

   call output_init( nvar, iprint_true )
   ! ^_^ TODO: move all cycle output to the output module; too painful right now

   write(xtime_fh,xtime_cycfmt) icountmod, agej, t9, rho, ONE-sum(yps(1,:)), ye, &
         (max(yps(1,iprint_true(i)),1.d-99), i = 1, nvar)
   call printonecycle(icountmod, dt_yr, agej, t9, rho, yps, acsp, csp, an, zn,.true.)

   ineut = ispe('NEUT ')

   if (print_rate) call print_one_rate(nvar, nvrel, yps, nu)

   ! cycle loop
   tsetup = wallclocktime() - tstart
   tstart = wallclocktime()
   modr: do while (.true.)

      select case( nsource )
      case(0)
         ! ^_^ hydrostatic (constant temperature and density)
         if ( agej >= t_max_yr ) then
            exit modr
         end if
         icountmod = icountmod + 1
         dt_yr     = min( dt_yr * dt_factor, dt_max_yr )
         dt        = dt_yr * yrs2sec
         agej      = agej + dt_yr
         t9_0 = t9 ; t9_1 = t9 ; rho_0 = rho ; rho_1 = rho
      case(1)
         ! trajectory
         if ( icountmod == t%ntimesteps - 1 ) then
            call printonecycle(icountmod, dt_yr, agej, t9, rho, yps, acsp, csp, an, zn, .true.)
            exit modr
         end if
         icountmod = icountmod + 1
         if ( t%timestep_given ) then
            dt    = t%time(icountmod)
            t9_0  = t%temp(icountmod); t9_1 = t%temp(icountmod)
            rho_0   = t%dens(icountmod); rho_1 = t%dens(icountmod)
         else
            dt        = t%time(icountmod + 1) - t%time(icountmod)
            t9_0      = t%temp(icountmod)
            t9_1      = t%temp(icountmod + 1)
            rho_0     = t%dens(icountmod)
            rho_1     = t%dens(icountmod + 1)
         end if
         ! arithmetic mean (instead of geometric mean, as was done previously...)
         t9           = HALF * (t9_0 + t9_1)
         rho          = HALF * (rho_0 + rho_1)
         dt_yr        = dt / yrs2sec
         agej         = agej + dt_yr
      end select

      call calculate_ye(yps(1,:), an, zn, ye, considerisotope)

      ! ^_^ physics
      call evaluate_all_rates(ye, nvar, nvrel, rho, t9, yps, nu)

      ! ^_^ solver
      call integrate_network( nvar, yps, t9_0, t9_1, rho_0, rho_1, ye, dt, nvrel, nu, ierr )
      
      select case( ierr )
      case( 0 )
         ! all fine
         continue
      case( MAX_SUBSTEPS )
         stop "max. substeps reached"
      case( SMALL_DT )
         stop "time step underflow"
      case default
         stop "exception raised in network integration"
      end select

      call printonecycle(icountmod, dt_yr, agej, t9, rho, yps, acsp, csp, an, zn, .false.)

      xndens = yps(1,ineut) * rho * avogadro

      ! STDOUT
      if (mod(icountmod,terminal_cycles) == 0) then
         write(*,stdout_fmt) icountmod, agej, xndens, t9, rho, ye, &
            ntime, ntime_last, time_matrix_inv, nvar1, iter, titers, nsubt
      end if

      ! *** print x-time.dat every time step
      write(xtime_fh, xtime_cycfmt) icountmod, agej, t9, rho, ONE-sum(yps(1,:)), ye, &
         (yps(1,iprint_true(i)), i=1,nvar)

   end do modr

   print *, "simulation complete .."
   write(*,"(a17,es10.3,a2)") "setup time       ", tsetup, " s"
   write(*,"(a17,es10.3,a2)") "computation time ", wallclocktime() - tstart, " s"

   call output_fin()

end program ppn

subroutine print_one_rate(nvar,nvrel,yps, nu)
   use nuc_data, only: zis
   use evaluate_rates, only: evaluate_all_rates
   use array_sizes
   use frame_knobs, only: which_rate
   use rates, only: v, k1, k2, k3, k4, k5, k6, k7, k8
   use utils
   use constants
   use neutrinos, only: neutrino
   implicit none
   integer(i4) :: nvar, nvrel
   real(r8) :: yps(nsp)
   integer(i4), parameter :: num = 1001
   real(r8), parameter :: &
      t9min = 0.01_r8, t9max = 10._r8, &
      rho = 1.e6_r8, ye = 0.5_r8
   integer(i4) :: i, fh, w
   real(r8) :: dlogt, t9, logt(num), rate(num)
   type(neutrino) :: nu

   ! temperature grid interval and boundary
   dlogt = (log10(t9max*1e9_r8) - log10(t9min*1e9_r8)) / (num - 1)
   logt(1) = log10(t9min*1e9_r8)

   ! open file for writing and write header
   print *, "writing rate for reaction ", which_rate, "to file"
   open( file = "rate.txt", newunit = fh, action = "write", status = "new")
   w = which_rate
   write(fh,"(2X,I2,2X,A5,'  + ',I2,2X,A5,2X,'->',1X,I2,2X,A5,'  + ',I2,2X,A5)") &
      k2(w), zis(k1(w)), k4(w), zis(k3(w)), k8(w), zis(k7(w)), k6(w), zis(k5(w))
   write(fh,*) "rho [g/cm^3] = ", rho 
   write(fh,"(a13,a23)") "T[GK]", "rate"

   do i = 1, num
      if (i > 1) logt(i) = logt(i-1) + dlogt
      t9 = 10._r8**logt(i)/1.e9_r8
      call evaluate_all_rates(ye, nvar, nvrel, rho, t9, yps, nu)
      write(fh,"(es13.6,es23.15)") t9, v(w)
   end do

   close(fh)
   print *, "rate for reaction ", w, "written to file rate.txt"
   stop

   end subroutine print_one_rate

!---
! DESCRIPTION:
!
!> @brief calculate the timescale of each reaction
!> @todo make sure that this is correct, and if it is, see whether it is useful
!> for a time step!estimator. Actually the calculate_rhs already returns this
!---
subroutine calculate_timescale(yps, yn, yp, ya, anetw, znetw, time_scale_reac)
      use array_sizes
      use nuc_data, only: considerreaction
      use rates
      use utils, only: r8
      use constants
      implicit none
      real(r8) :: yn, yp, ya, time_scale_reac(nre),yps(nsp)
      real(r8) :: anetw(nsp),znetw(nsp)
      integer i

      do i = 1, nre
         if (.not. considerreaction(i) .or. v(i) < 1.d-90) then
            time_scale_reac(i) = 1.e90_r8
            cycle
         end if
         if (int(anetw(k3(i))) == 1.and.int(znetw(k3(i))) == 1) then
            time_scale_reac(i) = ONE / (yp * v(i))
         else if (int(anetw(k3(i))) == 1.and.int(znetw(k3(i))) == 0) then
            time_scale_reac(i) = ONE / (yn * v(i))
         else if (int(anetw(k3(i))) == 4.and.int(znetw(k3(i))) == 2) then
            time_scale_reac(i) = ONE / (ya * v(i))
         else if (k2(i) == 1.and.k4(i) == 0) then
            time_scale_reac(i) = ONE / v(i)
         else if (k2(i) > 1.and.k4(i) == 0) then
            time_scale_reac(i) = ONE / ((yps(k1(i))/anetw(k1(i)))**(k2(i)-1) * v(i))
         else
            time_scale_reac(i) = ONE / (yps(k3(i))/anetw(k3(i))*v(i))
         end if
      end do
end subroutine calculate_timescale


!---
! DESCRIPTION:
!
!> @brief Writes iso_massf and flux output files
!---
subroutine printonecycle( icountmod, dzeit, agej, t9, rho, yps, acsp, csp, anetw, znetw, force_write )
      use array_sizes
      use nuc_data
      use rates
      use frame_knobs
      use utils, only: i4, r8, table_output
      use sorting, only: indexx
      use constants, only: avogadro
      use networksetup, only: nrnc1, nvnc1
      use reaction_info, only: ilabb, bind_energy_diff
      use physics_knobs, only: decay

      implicit none

      integer(i4) :: icountmod, i, modell, indxx(nre)
      real(r8), intent(in) :: dzeit, agej, t9, rho, yps(nsp), acsp(nsp)
      character(len=5) :: csp(nsp), cmodell
      real(r8) :: anetw(nsp),znetw(nsp)
      real(r8), dimension(nre) :: flux_to_print, time_scale_reac, flux_to_print_0
      logical :: force_write
      common /flux/ flux_to_print
      common /flux_0/ flux_to_print_0

      real(r8) :: energy_generated(nre)
      real(r8) :: total_energy_flux, yn, yp, ya

      total_energy_flux = 0._r8

      ! define yn, yp,ya and reaction timescales
      ! *** NOTE: rates are already multiplied by density term in the physics package.
      ! *** We do not need to multiply again here.
      yn = yps(ispe('NEUT '))
      yp = yps(ispe('PROT '))
      ya = yps(ispe('HE  4')) * 0.25_r8
      ! ***
      if (iplot_flux_option == 1 .and. icountmod > 0) then
         time_scale_reac = 0._r8
         call calculate_timescale(yps, yn, yp, ya, anetw, znetw, time_scale_reac)
      end if
      ! print each model in its own file
      ! checking icount
      modell = icountmod
      if (mod(modell,print_cycles)  ==  0 .or. decay .or. force_write) then
         write(cmodell,'(I5)')modell
         if (modell  <=  9 ) then
            cmodell="0000"//cmodell(5:5)
         else if (modell  <=  99 ) then
            cmodell="000"//cmodell(4:5)
         else if (modell  <=  999 ) then
            cmodell="00"//cmodell(3:5)
         else if (modell  <=  9999 ) then
            cmodell="0"//cmodell(2:5)
         endif
         if (decay) cmodell = "decay"
         ! *** printed nucleosynthesis fluxes, [dY/dt],
         ! *** calculated and printed and energy fluxes [erg/g]
         if (iplot_flux_option == 1 .and. icountmod > 0) then
            ! ***  using integrated fluxes or d_flux/dt?
            if (i_flux_integrated  ==  0) then
               flux_to_print_0 = flux_to_print
            else if (i_flux_integrated  ==  1) then
               flux_to_print_0 = flux_to_print_0 + flux_to_print
            end if
            open(newunit=table_output,file='flux_'//cmodell//".DAT")
            write(table_output,'(3x,a,2x,4(a,2x),6x,4(a,2x),2x,a,2(4x,a))')'#', &
                  'Z_k1','A_k1','Z_k3','A_k3','Z_k5','A_k5','Z_k7','A_k7', &
                  'flux [dY/dt]','energy [erg/(g*s)]','timescale reac [s]'
            do i = 1, nrnc1
               write(table_output,'(i5,4(2x,i3),8x,4(2x,i3),3(8x,ES12.5))') i, &
                     int(znetw(k1(i))),int(anetw(k1(i))), &
                     int(znetw(k3(i))),int(anetw(k3(i))), &
                     int(znetw(k5(i))),int(anetw(k5(i))), &
                     int(znetw(k7(i))),int(anetw(k7(i))), &
                     max(1.0d-99,flux_to_print_0(i)), &
                     max(1.0d-99,flux_to_print_0(i)*bind_energy_diff(i)), &
                     max(1.0d+99,time_scale_reac(i))
            end do

            ! here below the total energy flux is calculated
            !             total_energy_flux = 0.
            energy_generated(1:nre) = flux_to_print_0(1:nre) * bind_energy_diff(1:nre)
            total_energy_flux = sum(energy_generated)

            if (iolevel  >  2) then
               ! *** if I want to know the 20 highest energy fluxes
               ! *** I do this, with high IO
               write(*,*) ' '
               write(*,*) 'top 20 energy fluxes:'
               call indexx(nre,energy_generated,indxx)
               write(*,03) (ZIS(K1(indxx(i))),ZIS(K3(indxx(i))), &
                     energy_generated(indxx(i)), i=nre,nre-19,-1)
               write(*,*) ' '
            end if

            close(table_output)        
         end if

         ! *** iso_massf is written
         open(newunit = table_output,file=cprefix(1:len_trim(cprefix))//cmodell//".DAT")
         write(table_output,'(A1,1x,A3,2x,2(3X,A1,1x),A4,2x,A12,2x ,A5)') &
               'H','NUM','Z','A','ISOM','ABUNDANCE_MF','ISOTP'
         write(table_output,'(1x,A,1x,I8,1x,A,1x,ES9.2,1x,A,1x,ES11.4)') &
               '# mod',icountmod,'dzeit',dzeit,'agej',agej
         write(table_output,'(2(1x,A,2x,ES10.3))') '# t9 = ', t9, 'rho = ', rho
         write(table_output,'(1x,A,1x,ES12.5)') '# densn', yn * rho * avogadro
         write(table_output,'(1x,A,1x,ES12.5)') '# densp', yp * rho * avogadro
         write(table_output,'(1x,A,1x,ES12.5)') '# densa', ya * rho * avogadro
         write(table_output,'(1x,A,1x,ES12.5)') '# total_energy_flux_[erg/(g*s)]', total_energy_flux

         do i=1,nsp
            if(considerisotope(i))then
               write(table_output,'(I5,2x,2(1X,F4.0),2x,I2,2x,ES12.5,2x,A5)') i, &
                     znetw(i), anetw(i), isomeric_state(i), max(yps(i),1.0d-99), zis(i)
            end if
         end do
         close(table_output)    
      end if

      03    format(1x,3(2(a5),'=',1pe10.3,'  '))

end subroutine printonecycle


!---
! DESCRIPTION:
!> @brief Prints splash screen
!---
subroutine print_splash
      use utils, only: get_unused_file_handle
      integer :: fh
      fh = get_unused_file_handle()

      write(*,*) ''
      write(*,*) ''
      write(*,*) '       ooooooooo.   ooooooooo.   ooooo      ooo '
      write(*,*) '       `888   `Y88. `888   `Y88. `888b.     `8` '
      write(*,*) '        888   .d88`  888   .d88`  8 `88b.    8  '
      write(*,*) '        888ooo88P`   888ooo88P`   8   `88b.  8  '
      write(*,*) '        888          888          8     `88b.8  '
      write(*,*) '        888          888          8       `888  '
      write(*,*) '       o888o        o888o        o8o        `8  '
      write(*,*) ''
      write(*,*) ''
      write(*,*) ''
end subroutine print_splash


