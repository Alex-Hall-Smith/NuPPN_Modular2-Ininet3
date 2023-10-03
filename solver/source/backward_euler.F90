module backward_euler
   use frame_knobs
   use rates
   use array_sizes
   use utils
   use jac_rhs
   use nuc_data
   use solver_knobs
   use constants
   use sparsekit_module, only: dump_mtx, coocsr
   use linalg, only: solve_Ax_equals_b_csr
   use solver_diagnostics
   implicit none
   private

   public :: backward_euler_newton, newton_iteration

contains

   !---
   ! DESCRIPTION:
   !
   !> @author Falk Herwig, Marco Pignatari, Samuel Jones
   !> @brief This subroutine performs a Newton-Raphson iteration scheme for a fully-implicit, first-order accurate network integration
   ! (backward Euler).
   ! 
   !> @details Derived from: networx99.f,v 1.7 1999/05/20 09:47:14 fherwig Exp fherwig (c) Falk Herwig, 1999, 2007
   !---
   subroutine backward_euler_newton( nvar, yps, t9, rho, dt, dtrcmd, ierr )
         use errors
         integer :: &
               nvar, &       !< number of isotopes in full network
               igrd, &       !< index of isotope with greatest relative correction after integration step
               ierr
         real(r8) :: &
               yps(nsp), &     !< isotopic mass fractions
               dt, &           !< network time step (seconds)
               y_old(nsp), &   !< converged abundances at previous time step
               y_cur(nsp), &   !< abundances at current iteration
               y_new(nsp), &   !< abundances calculated during new iteration
               grd, &          !< greatest relative difference from previous iteration
               dtrcmd, &       !< recommended time step
               err, &          !< integration error
               tol, &          !< integration error tolerance
               rho, t9, sum_check, tti1, tti2
         real(r8), parameter :: &
               rmin = 0.2_r8, rmax = 2

         titers = 0 ; ntime = ZERO ; nvar1 = 0 ; iter = 0
         sum_check  = grdthreshold * nvar * dgrd

         y_new(:) = yps(:); y_old(:) = yps(:); y_cur(:) = yps(:)
         ierr = 0

         newton: do                        
            ! if convergence is not reached in ittd iterations, recommend smaller time step
            ! and exit
            if ( iter > ittd ) then 
               if (maxval(-y_new) >= sum_check) then
                  ! solution converged with "large" negative abundances, error is no longer grd
                  ! but the amount by which the min. abundance is below the -SMALL threshold
                  err = maxval(-y_new); tol = sum_check
                  ierr = NEG_ABUNDS
               else if (abs(ONE - sum(y_new)) > SUM_TOL ) then
                  ! mass conversation poor
                  err = abs(ONE - sum(y_new)); tol = SUM_TOL
                  ierr = SUM_X
               else
                  ! newton-raphson did not converge to requested accuracy
                  err = grd; tol = dgrd
                  ierr = SLV_CNVRG
               end if
               ! recommended time step to try instead
               dtrcmd = max(min(( tol / max(err,SMALL) ) ** HALF, rmax), rmin) * dt
               return
            end if

            ! perform one newton-raphson iteration
            tti1 = wallclocktime()
            call newton_iteration( y_old, y_cur, dt, grd, igrd, nvar, y_new )
            ntime = wallclocktime() - tti1
            if (iolevel >= 3) write (*,fmtprofile) ' NRNW cpu_time: ', ntime

            iter = iter + 1 ; titers = titers + 1

            ! Convergenge check (also on the negative values; possible solution)
            ! If a negative value is higher than sum_check, step still fails
            ! If mass conservation is violated, step still fails
            if ( grd <= dgrd .and. maxval(-y_new) < sum_check .and. abs(ONE - sum(y_new)) <= SUM_TOL ) then
               ! adequately converged
               yps(:) = y_new(:)

               if ( iolevel >= 3 ) then
                  print *, "successful step"
                  write(*,'(A12,I3)') "iterations= ", iter
                  write(*,'(3A,ES12.4,A,ES12.4)') "Isotope ", zis(igrd), " gr.rel.Cor=", grd, " with abund.=", y_new(igrd)
               end if

               ! recommended time step to try next
               dtrcmd = TWO*max(min(( dgrd / max(grd,SMALL) ) ** HALF, rmax), rmin) * dt
               call ypsplus(yps)
               exit newton
            end if    

            y_cur(:) = y_new(:)
         end do newton

         ntime = ntime / iter
   end subroutine backward_euler_newton


   !---
   ! DESCRIPTION:
   !
   !> @brief This routine performs one Newton-Raphson iteration, called from nucnet99
   !
   !> @details Id: nrnw.f,v 1.3 1999/05/20 09:47:27 fherwig Exp fherwig
   !> (c) Falk Herwig, 1999, 2007; Samuel Jones 2015-2017
   !---
   subroutine newton_iteration( y_old, y_cur, dt, grd, igrd, nvar, y_new )
         real(r8), dimension(nsp) :: &
               y_old, &    !< isotopic mass fractions at end of previous time step
               y_cur, &    !< isotopic mass fractions at last iteration
               y_new, &
               delta_y, &  !< correction to solution for previous iteration
               rhs, &      !< the network right hand side
               grdarr
         real(r8) :: dt, grd, tti1, tti2, dt_est
         integer :: nvar, igrd, i, ice
         logical :: mask(nsp)

         tti1 = wallclocktime()
         call calculate_jacobian( dt, mask, y_cur )
         time_jacobian = wallclocktime() - tti1

         jac_coo % val = jac_coo % val * ( -dt )
         call reduce_jacobian( mask, nvar1 )

         ! calculate RHS, compress for problem size
         call calculate_dxdt( y_cur, dt_est )

         ! compress for problem size
         rhs = y_old(:) - y_cur(:) + dt * dxdt(:)
         call reduce_rhs( nvar1, mask, y_new, y_cur, rhs )

         !> solve the system
         tti2 = wallclocktime()
         call solve_Ax_equals_b_csr( reduced_jac_csr, reduced_rhs )
         tti1 = wallclocktime()
         time_matrix_inv = tti1 - tti2

         ! apply correction. 
         reduced_y_new(1:nvar1) = reduced_y_cur(:) + reduced_rhs(:)

         ! write arrays back to length nsp (reduced_y_new ---> y_new).
         ! y_new is the new abundance calculated in the present interation.
         call ddecompac( mask, reduced_y_new, y_new, 1, nsp, nvar1, y_new )

         delta_y(:) = y_new(:) - y_cur(:)

         ! GRD == "greatest relative difference" from previous iteration
         grdarr(:) = ZERO
         where ( y_cur(:) > grdthreshold )
            grdarr(:) = abs( delta_y(:) / y_cur(:) )
         end where
         grd = maxval( grdarr ); igrd = maxloc( grdarr , dim = 1)

         if (iolevel >= 3) then
            write( *, fmtprofile ) "nrnw_jacobian_time/s:"     ,     time_jacobian
            write( *, fmtprofile ) "nrnw_mtx_inv_time/s:"      ,     time_matrix_inv
            write(*,'(3(1X,A),1X,ES9.2,1X,A,1X,ES9.2)') "Rel. GKOR", zis(igrd), "GRD=", grd, "X(JGRD)=", y_cur(igrd)
         end if

         call free_reduced_jacobian()
         call free_reduced_rhs()
   end subroutine newton_iteration


end module backward_euler
