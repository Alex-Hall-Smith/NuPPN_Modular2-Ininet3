module solver_knobs
   use array_sizes
   use utils
   real(r8) :: &
         grdthreshold, &    !< minimum abundance for checking relative correction to check for convergence
         cyminformal, &     !< the solver zero
         dgrd, &            !< convergence limit for largest relative correction
         min_t9_for_slv, &  !< minimum temperature in GK for solving the network
         delta_state_max, & !< max. allowed relative change in state (T, rho, ye) over a time step
         min_rate_ov_dt     !< lower limit for reaction as fraction of time step in order to use it in dynamic network

   real(r8), parameter :: &
         ye_tol = 1.e-5_r8, &  !< (relative) error tolerance for coupling of weak reactions in NSE
         SUM_TOL = 1.e-6_r8  !< mass conservation tolerance (deviation from 1)

   integer(i4), parameter :: &
         EULER  = 0, &
         DEUFL  = 1, &
         RK4    = 1

   integer :: &
         ittd, &               !< max N-R iterations
         irdn, &               !< switch for static (0?) or dynamic (2) network
         mat_solv_option, &    !< choice of matrix inversion technique
         nse_max_wk_iters, &   !< max. allowed iterations in weak NSE weak coupling
         max_sub_steps, &      !< max. allowed sub time steps
         integration_method, & !< 0 (EULER - backward euler) or 1 (DEUFL - bader-deuflhard)
         nse_weak_method       !< 0 (EULER - forward euler) or 1 (RK$ - cash-karp rk4)

   logical :: &
         interpolate_temp, & !< interpolate temperature during sub time steps
         interpolate_dens    !< interpolate density during sub time steps

contains
   !> @brief Read ppn_solver.input file
   subroutine readsolverinput()
         use communication
         integer :: fh
         namelist/ppn_solver/ ittd, dgrd, grdthreshold, irdn, cyminformal, mat_solv_option, &
               min_rate_ov_dt, nse_max_wk_iters, interpolate_temp, interpolate_dens, integration_method, &
               max_sub_steps, nse_weak_method, min_t9_for_slv, delta_state_max

         ittd               = 6
         dgrd               = 1.e-3_r8
         grdthreshold       = 1.e-15_r8
         irdn               = 2
         cyminformal        = 1.0e-30_r8
         min_rate_ov_dt     = 1.e-99_r8
         mat_solv_option    = 5
         nse_max_wk_iters   = 1000
         interpolate_temp   = .true.
         interpolate_dens   = .true.
         integration_method = EULER
         nse_weak_method    = RK4
         max_sub_steps      = 4000
         min_t9_for_slv     = 0.01_r8
         delta_state_max    = 1.e99_r8

         if (master) then
            fh = get_unused_file_handle()
            open( fh, file = "ppn_solver.input" )
            read( fh, NML = ppn_solver )
            close(fh)
         end if

         if (cyminformal < 1.e-99_r8) then
            if (master) print *, 'full, static network in use: irdn choice irrelevant.'
            min_rate_ov_dt = -1.e99_r8
         end if

         if (irdn <= 0 .or. irdn > 3) stop 'invalid irdn choice'

#ifndef PPN
         call broadcast(ittd)             ; call broadcast(dgrd)
         call broadcast(irdn)             ; call broadcast(cyminformal)
         call broadcast(mat_solv_option)  ; call broadcast(nse_max_wk_iters)
         call broadcast(interpolate_dens) ; call broadcast(integration_method)
         call broadcast(max_sub_steps)    ; call broadcast(grdthreshold)
         call broadcast(min_rate_ov_dt)   ; call broadcast(interpolate_temp)
         call broadcast(nse_weak_method)  ; call broadcast(delta_state_max)
#endif
   end subroutine readsolverinput


end module solver_knobs
