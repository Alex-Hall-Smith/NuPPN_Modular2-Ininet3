! MODULE: solver
!
!> @author Samuel Jones, Falk Herwig, Marco Pignatari
!
!> @todo RDN functionality should not go into jacobn

#undef VERBOSE

module solver
   use array_sizes, only: nsp, nre
   use utils, only: r8, i4
   use constants, only: ZERO, HALF, SMALL
   use jac_rhs
   use solver_knobs

   implicit none
   private

   logical :: do_nse = .false.

   public solver_init, integrate_network
#ifndef PPN
   public solver_broadcasts
#endif

contains

   !> @brief Initialization of the solver module
   subroutine solver_init()
         call jac_rhs_init()
   end subroutine solver_init


#ifndef PPN
   !> @brief Broadcasting some solver variables to the slaves
   subroutine solver_broadcasts()
         use communication
         return

         ! ^_^ integers
         call broadcast(mat_solv_option)
         call broadcast(irdn)
         call broadcast(ittd)

         ! ^_^ dp reals
         call broadcast(grdthreshold)
         call broadcast(dgrd)
         call broadcast(cyminformal)
   end subroutine solver_broadcasts
#endif

   !> todo add doxygen comments
   subroutine integrate_network( nvar, yps, t9_0, t9_1, rho_0, rho_1, ye, dtwant, nvrel, nu, ierr )
         use interp, only: create_1d_interpolant, interpolate_1d, interpolant
         use physics_knobs, only: i_nse, nse_min_t9
         use solver_diagnostics, only: nsubt, ntime, iter
         use errors, only: MAX_SUBSTEPS, SMALL_DT, error_code
         use nse_solver, only: nse_solver_cache_lm_ov_a
         use evaluate_rates, only: evaluate_all_rates
         use jac_rhs, only: calculate_dxdt
         use nuc_data, only: ispe, zis
         use neutrinos, only: neutrino
         use rates, only: v
         use communication
         integer :: nvar, nvrel, ierr, i, j, nsucc, nfail, err_fh
         real(r8) :: &
               yps(nsp), yps_old(nsp), yps_ent(nsp), ye, ye_old, &
               dtwant, dttry, dtdone, dtrcmd, dt_est, &
               t_0, t_1, &
               t9_0, t9_1, my_t9_0, my_t9_1, my_t9_last, my_t9, &
               rho_0, rho_1, my_rho_0, my_rho_1, my_rho_last, my_rho, &
               nw, thn, tti1, diff
         type(interpolant) :: t9intrp, rhointrp
         type(neutrino) :: nu
         logical :: frst = .true.
         character(len=5) :: lab(nre), labb(nre)
         common/lab/lab,labb

         dtdone = ZERO ; dtrcmd = ZERO ; yps_old(:) = yps(:) ; ye_old = ye
         yps_ent(:) = yps(:)
         my_t9_last = HALF * ( t9_0 + t9_1 ) ; my_rho_last = HALF * ( rho_0 + rho_1 )
         nsucc = 0 ; nfail = 0 ; iter = 0

         if ( interpolate_temp ) then
            ! set up temperature interpolant over step from time = 0 to dtwant
            call create_1d_interpolant ( 2 , [ZERO , dtwant] , [t9_0  , t9_1]  , t9intrp  ) 
         end if
         if ( interpolate_dens ) then
            ! set up density interpolant over step from time = 0 to dtwant
            call create_1d_interpolant ( 2 , [ZERO , dtwant] , [rho_0 , rho_1] , rhointrp ) 
         end if


         if ( .false. ) then
            ! estimate appropriate time step from max. flux
            call calculate_dxdt( yps, dt_est )
            dttry = min( dt_est, dtwant )
         else
            ! start by trying time step that was requested
            dttry = dtwant
         end if

         ! profiling
         tti1 = wallclocktime()


         ! sub time step loop
         time_step: do i = 1, max_sub_steps
            if ( i > 1 ) frst = .false.
            
            if ( interpolate_temp ) then
               ! get the average temperature over the current interval using interpolation
               t_0 = dtdone ; t_1 = min( dtdone + dttry, dtwant )
               call interpolate_1d ( t9intrp  , t_0 , my_t9_0  ) 
               call interpolate_1d ( t9intrp  , t_1 , my_t9_1  ) 

               ! use average temperature for the time step
               my_t9  = HALF * ( my_t9_0  + my_t9_1  ) 
            else
               ! use average temperature for the time step
               my_t9  = HALF * ( t9_0  + t9_1  ) 
            end if

            if ( interpolate_dens ) then
               ! get the average density over the current interval using interpolation
               t_0 = dtdone ; t_1 = min( dtdone + dttry, dtwant )
               call interpolate_1d ( rhointrp , t_0 , my_rho_0 ) 
               call interpolate_1d ( rhointrp , t_1 , my_rho_1 ) 

               ! use average density for the time step
               my_rho = HALF * ( my_rho_0 + my_rho_1 ) 
            else
               ! use average density for the time step
               my_rho = HALF * ( rho_0 + rho_1 ) 
            end if

            ! if temperature or density change by more than delta_state_max
            ! (relative change) -> smaller sub time step
            diff = abs(my_t9_0 - my_t9_1)/my_t9_0
            diff = max(diff, abs(my_rho_0 - my_rho_1)/my_rho_0)
            if (diff > delta_state_max .and. t9_0 > min_t9_for_slv .and. t9_1 > min_t9_for_slv) then
               dttry = HALF*dttry
               cycle
            end if

            if ( my_t9 /= my_t9_last .or. my_rho /= my_rho_last ) then
               ! temperature or density (or both) changed, so update rates
               call evaluate_all_rates( ye, nvar, nvrel, my_rho, my_t9, yps, nu )
               my_rho_last = my_rho ; my_t9_last = my_t9
            end if

            ! check if we should do NSE
            if ( i_nse == 1 .and. my_t9 >= nse_min_t9 ) then
               do_nse = .true.
               call nse_solver_cache_lm_ov_a()
            else
               do_nse = .false.
            end if

            ierr = 0

            if ( my_t9 >= min_t9_for_slv ) then
               ! perform time step
               call do_one_time_step( nvar, yps, my_t9, my_rho, ye, dttry, dtrcmd, frst, ierr )
               if ( dttry < SMALL ) then
                  print *, "process ", ipid, "failed"
                  write(*,"(a10,es23.15)") "t9_0 = ", t9_0
                  write(*,"(a10,es23.15)") "t9_1", t9_1
                  write(*,"(a10,es23.15)") "rho_0 = ", rho_0
                  write(*,"(a10,es23.15)") "rho_1 = ", rho_1
                  write(*,"(a10,es23.15)") "dt = ", dtwant
                  write(*,"(a10,es23.15)") "dtrcmd = ", dtrcmd
                  write(*,"(a10,es23.15)") "dttry = ", dttry
                  call write_problem(dtwant, t9_0, t9_1, rho_0, rho_1, ye, yps_ent)
                  print *, "SMALL_DT: check problem.txt"

                  ierr = SMALL_DT
                  return
               end if
            else
               ! no nucleosynthesis assumed, pass over time step
               ierr = 0
            end if

            ! check success of time step
            if ( ierr == 0 ) then
               ! step was successful, add the time done and check for completion of step
               dtdone = dtdone + dttry
               yps = max( yps, 1.e-99_r8 )
               yps_old(:) = yps(:)
               ye_old = ye
               nsucc = nsucc + 1
               if ( dtdone >= dtwant ) exit time_step
            else
               ! step failed, reset abundances
               yps(:) = yps_old(:)
               ye = ye_old
               nfail = nfail + 1
            end if

#if defined (VERBOSE) || defined(PPN)
            if ( mod( i, 20 ) == 0 ) then
               write(*,"(a14,i6,a1,3x,a8,es9.2,3x,a8,es9.2,3x,a8,es9.2,1x,a3,i5,1x,a3,i5,1x,a3,i3)") &
                     "sub time step ", i, ':', &
                     "dt wnt: ", dtwant, &
                     "dt act: ", dttry, &
                     "dt dne: ", dtdone, &
                     "s: ", nsucc, &
                     "f: ", nfail, &
                     "e: ", ierr
            end if
#endif

            ! new time step to try
            dttry = min( dtrcmd, dtwant - dtdone )

            ! check for max. substeps reached, and return error
            if ( i == max_sub_steps ) then
               print *, "process ", ipid, "failed; t9_0 = ", t9_0, "t9_1", &
                  t9_1, "rho_0 = ", rho_0, "rho_1 = ", rho_1, "dt = ", dtwant
               print *, "MAX_SUB_STEPS: check problem.txt"
               call write_problem(dtwant, t9_0, t9_1, rho_0, rho_1, ye, yps_ent)
               ierr = MAX_SUBSTEPS
               return 
            end if
         end do time_step


         ! diagnostics
         nsubt = i ; ntime = wallclocktime() - tti1

         ! if using Bader-deuflhard solver, mass conservation should be preserved to 10% of SUM_TOL,
         ! so renormalise here
         if (integration_method == DEUFL) yps(:) = yps(:) / sum(yps)

         ! clip abundances
         yps(:) = max( yps(:), 1.e-99_r8 )
   end subroutine integrate_network


   subroutine do_one_time_step( nvar, yps, t9, rho, ye, dt, dtrcmd, frst, ierr )
         use physics_knobs, only: i_nse, nse_min_t9
         use nuc_data, only: considerreaction, considerisotope
         use array_sizes, only: nre
         use nse_solver, only: nse_yedot_rk4, nse_yedot_euler
         use bader_deuflhard, only: bader_deuflhard_integration
         use backward_euler, only: backward_euler_newton
         integer :: nvar
         integer, intent(out) :: ierr
         real(r8), intent(in) :: yps(nsp), t9, dt, rho, ye
         real(r8), intent(out) :: dtrcmd !< recommended time step from solve
         real(r8) :: an(nsp), zn(nsp)
         logical :: frst
         common/ cnetw / an, zn

         ierr = 0

         if ( do_nse ) then
            select case ( nse_weak_method )
            case( EULER )
               call nse_yedot_euler(nvar, yps, t9, rho, ye, dt, dtrcmd, ierr)
            case( RK4 )
               call nse_yedot_rk4(nvar, yps, t9, rho, ye, dt, dtrcmd, ierr )
            end select
         else
            select case ( integration_method )
            case ( EULER )
               call backward_euler_newton( nvar, yps, t9, rho, dt, dtrcmd, ierr )
            case ( DEUFL )
               call bader_deuflhard_integration( yps, dt, dtrcmd, frst, ierr )
            end select
            call calculate_ye( yps, an, zn, ye, considerisotope )
         end if
   end subroutine do_one_time_step

   !> write information about a time step that failed
   subroutine write_problem(dt, t9_0, t9_1, rho_0, rho_1, ye, yps_ent)
      use nuc_data, only: zis, considerisotope, considerreaction
      use rates, only: v
      use reaction_info, only: lab
      implicit none
      integer(i4) :: i, fh
      real(r8) :: dt, t9_0, t9_1, rho_0, rho_1, ye, yps_ent(nsp)

      !> write conditions
      open(file="problem.txt", newunit=fh)
      write(fh,"(a10,es23.15)") "t9_0", t9_0
      write(fh,"(a10,es23.15)") "t9_1", t9_1
      write(fh,"(a10,es23.15)") "rho_0", rho_0
      write(fh,"(a10,es23.15)") "rho_1", rho_1
      write(fh,"(a10,es23.15)") "ye", ye
      write(fh,"(a10,es23.15)") "dt / s", dt

      write(fh,*)

      !> write abundances
      do i = 1, nsp
         if (.not. considerisotope(i)) cycle
         write(fh,"(a5,es23.15)") zis(i), yps_ent(i)
      end do

      !> write reaction rates
      do i = 1, nre
         if (.not. considerreaction(i)) cycle
         write(fh,"(i9,es23.15,a8)") i, v(i), lab(i)
      end do

      close(fh)

   end subroutine write_problem

end module solver
