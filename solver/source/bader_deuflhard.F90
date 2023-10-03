#undef VERBOSE

!> @author Samuel Jones

!> This module is an implementation of the Bader-Deuflhard time integration algorithm for a
!> system of differential equations. The method is a generalisation of the Bulirsch-Stoer
!> algorithm. It is a semi-implicit midpoint method of variable order that is combined with
!> Richardson extrapolation.

module bader_deuflhard
   use frame_knobs
   use array_sizes
   use utils
   use jac_rhs
   use nuc_data
   use solver_knobs
   use constants
   use sparsekit_module, only: dump_mtx, coocsr, amux
   use linalg, only: solve_Ax_equals_b_csr
   use solver_diagnostics
   use errors
   use communication
   implicit none
   private

   integer, parameter :: &
         FIRST = 1, MID = 2, LAST = 3, nlevels = 8

   integer(i4) :: &
         nseq(nlevels+1), & !< number of steps in each level of Bader-Deuflhard sequence
         kmax, &            !< max. row number used in the Richardson extrapolation
         kopt               !< optimal column for convergence

   real(r8) :: &
        qcol(nsp,nlevels), &      !< diagonal of extrapolation tableau
        yscal(nsp), &             !< error scaling
        yerr(nsp,nlevels), &      !< error
        err(nlevels), &           !< extrapolation error
        a(nlevels+1), &           !< work estimate for reaching T_ii in extrapolation tableau
        alf(nlevels,nlevels), &   !< correction factors
        errmax, &                 !< max. error of all equations
        abuerr, &                 !< error in abundance equation solution
        negerr, &                 !< error due to negative abundances
        sumerr, &                 !< mass conservation residual
        xest, x(nlevels)

  real(r8), parameter :: &
        cf  = 1, &            !< weight work needed for rhs eval.
        cj  = 0._r8 * cf, &   !< weight work needed for jacob eval.
        cs  = 10, &           !< weight work needed for forward/backward substitution triangular solve
        clr = cj, &           !< weight work needed for LU decomposition
        safe1 = FOURTH, &     !< safety factor
        safe2 = 0.7_r8, &     !< safety factor
        redmin = 0.7_r8, &    !< min. amount to modify time step by
        redmax = 1.e-5_r8, &  !< max. amount to modify time step by
        scalmx = 0.25_r8      !< 1/scalmx is the max. factor the time step can be increased by

  logical :: &
        use_weights = .true. !< use work weights to calculate work estimates

  public :: bader_deuflhard_init, bader_deuflhard_integration

contains

   !> initialisation of variables that are to be held constant over the simulation
   subroutine bader_deuflhard_init
         integer :: i
         real(r8) :: alpha = 5._r8 / 7._r8

         nseq(1:3) = [ 2, 6, 10 ]
         do i = 4, nlevels+1
            nseq(i) = nint( real(nseq(i-1),r8) / alpha )
         end do

         call compute_work_estimates
         call compute_correction_factors
         call determine_optimal_convergence
   end subroutine bader_deuflhard_init


   !> compute estimates of the computational work associated with reaching each column of he
   !> extrapolation tableau
   subroutine compute_work_estimates
         integer :: i

         if ( use_weights ) then
            a(1) = cj + clr + (nseq(1) + 1) * (cf + cs)
            do i = 2, nlevels
               a(i) = a(i-1) + nseq(i) * (cf + cs) + clr + cs
            end do
         else
            a(1) = nseq(1) + ONE
            do i = 2, nlevels + 1
               a(i) = a(i-1) + nseq(i)
            end do
         end if
   end subroutine compute_work_estimates


   !> compute the correction factors \f$ \alpha\ f$
   subroutine compute_correction_factors
         integer :: k, iq

         do iq = 2, nlevels
            do k = 1, iq - 1
               alf(k,iq) = ( dgrd * safe1 ) ** ( ( a(k+1) - a(iq+1) ) / ( ( a(iq+1) - a(1) + ONE ) * (TWO*k + 1) ) )
            end do
         end do
   end subroutine compute_correction_factors


   !> determine the row number kopt in which it will be most computationally efficient to
   !> achieve convergence in, considering computational effort of increasing the order but
   !> doing a larger time step vs doing more time steps but converging at lower order
   subroutine determine_optimal_convergence
         integer :: i

         do kopt = 2, nlevels - 1
            if ( a(kopt + 1) > a(kopt)*alf(kopt-1,kopt) ) then
               kmax = kopt
               exit
            end if
         end do
   end subroutine determine_optimal_convergence



   !> @brief This subroutine attempts to perform a semi-implicit extrapolation time integration
   !> of the reaction network over the time step dt ("Bader-Deuflhard" method).
   !> Should the step fail for reasons that could include large negative abundances,
   !> violation of mass conservation, non-convergence, an error flag will be returned
   !> indicating the issue along with a suggestion for an appropriate time step. If the
   !> solution converges, the algorithm will suggest a new, larger time step if increasing
   !> the time step would be computationally efficient
   !
   !> @todo interpolate temperature in time during substeps. This will probably require
   !> the use of high-order numerical derivatives because of the computational overhead
   !> of evaluating the reaction rates.
   subroutine bader_deuflhard_integration( yps, dt, dtrcmd, frst, ierr )
         use solver_diagnostics, only: iter
         integer :: &
               k, & !< row number in the extrapolation tableau
               kk, km, i, ierr
         real(r8) :: &
               yps(nsp), &   !< isotopic mass fractions
               dt, &         !< network time step (seconds)
               y_k(nsp), &   !< abundances at current sub step
               work, workmin, scal, fact, red, y_extrap(nsp), dtrcmd
         logical :: &
               mask(nsp), &     !< .true. where isotope should be included based on fluxes
               frst, &          !< is this a completely new integration (.true.) or just a sub step of the previous one?
               reduce = .false. !< whether or not to reduce the time step

         nvar1 = 0 ; scal = ONE ; km = -1 ; red = -ONE ; y_k(:) = yps(:); ierr = 0

         call calculate_jacobian( dt, mask, yps )
         jac_cpy%val(:) = jac_coo%val(:)

         do k = 1, kmax
            ! deuflhard sequence loop

            y_k(:) = yps(:)

            ! compute integration over number of sub steps in Bader-Deuflhard sequence level k
            call semi_implicit_extrapolation( dt, k, y_k, mask, ierr )

            ! check for terrible solutions and recommend reducing time step by a factor of 5
            if (ierr /= 0) then
               dtrcmd = dt * FIFTH
               return
            end if

            iter = iter + nseq(k)

            ! add to extrapolation tableau and attempt polynomial extrapolation to solution
            call extrapolate( k, y_k, y_extrap )

            ! update current solution with Richardson-extrapolated abundances
            y_k(:) = y_extrap(:)


            ! check for terrible solutions and recommend reducing time step by a factor of 5
            if ( any( y_k(:) /= y_k(:) ) ) then
               ierr = GOT_NAN
               dtrcmd = dt * FIFTH
               return
            end if

            ! error scaling - clip to grdthreshold instead of cyminformal
            ! (1e-15 instead of 1e-30 typically)
            !yscal(:) = max(abs(y_k(:)), cyminformal)
            yscal(:) = max(abs(y_k(:)), grdthreshold)

            ! compute and save max. error for this order
            if ( k /= 1 ) then
               call compute_error( k, y_k )
               km = k - 1
               ! max. error in extrapolated abundances y_k
               err(km) = ( errmax / safe1 ) ** ( ONE / (TWO * km + ONE ) )
            end if

            ! check for convergence if we are at level 2 or more and have either reached the optimal column
            ! for convergence in the extrapolation tableau for convergence or we are trying a new integration
            if ( k /= 1 .and. ( k >= kopt - 1 .or. frst ) ) then

               if ( errmax <= ONE ) then
                  ! converged: flag to try increasing time step and exit
                  reduce = .false.
                  exit
               else
                 ! not converged: check the circumstance and flag to reduce time step
                 if ( k == kmax .or. k == kopt + 1 ) then
                    ! reached max. order or overshot optimal convergence column without converging
                    reduce = .true. ; red = safe2 / err(km)

                 else if ( k == kopt ) then
                    ! reached optimal convergence column without converging
                    if ( alf(kopt-1,kopt) < err(km) ) then
                       reduce = .true. ; red = safe2 / err(km)
                    end if

                 else if ( kopt == kmax ) then
                    ! optimal convergence column is set at max. allowed order/level
                    if ( alf(km,kmax-1) < err(km)) then
                       reduce = .true. ; red = alf(km,kmax-1) * safe2 / err(km)
                    end if

                 else if ( alf(km,kopt) < err(km) ) then
                    reduce = .true. ; red = alf(km,kopt-1) / err(km)
                 end if
              end if

            end if
         end do

         if ( reduce ) then
            red = max(redmax, min(red, redmin))
            dtrcmd = dt * red
            !> set correct error flag
            if (negerr > sumerr .and. negerr > abuerr) then
               ierr = NEG_ABUNDS
            else if (sumerr > abuerr) then
               ierr = SUM_X
            else
               ierr = SLV_CNVRG
            end if
         else
            ! compute optimal row for convergence
            workmin = 1.e35_r8
            do kk = 1, km
               fact = max(err(kk),scalmx)
               work = fact * a(kk+1)
               if ( work < workmin ) then
                  scal    = fact
                  workmin = work
                  kopt    = kk + 1
               end if
            end do
            dtrcmd = TWO * dt / scal ! add factor two

            ! check for possible order increase, but not if step-size was just reduced.
            if ( kopt >= k .and. kopt /= kmax ) then
               fact = max(scal/alf(kopt-1,kopt),scalmx)
               if ( a(kopt+1) * fact <= workmin ) then
                  dtrcmd = dtrcmd * scal / fact
                  kopt = kopt + 1
               endif
            endif
         end if

         yps(:) = y_k(:)
   end subroutine bader_deuflhard_integration



   subroutine semi_implicit_extrapolation( time_step, k, y_k, mask, ierr )
         real(r8) :: time_step, h, y_k(nsp), delta_k(nsp)
         integer(i4), intent(in) :: k
         integer(i4) :: i, ierr
         logical :: mask(nsp)
         ierr = 0
         
         h = time_step / real(nseq(k),r8)
         xest = ( h / nseq(k) ) ** 2

         ! apply prefactor of substep h to jacobian and reduce (1-hJ) to problem size
         jac_coo%val(:) = jac_cpy%val(:) * (-h)
         call reduce_jacobian( mask, nvar1 )
         call convert_reduced_jac_to_csc()

         ! first step
         call bader_deuflhard_step( y_k, delta_k, h, nvar1, mask, FIRST, ierr )
         if (ierr /= 0) then
            call free_jac
            return
         end if

         ! mid steps
         do i = 1, nseq(k) - 1
            call bader_deuflhard_step( y_k, delta_k, h, nvar1, mask, MID, ierr )
            if (ierr /= 0) then
               call free_jac
               return
            end if
         end do

         ! end step
         call bader_deuflhard_step( y_k, delta_k, h, nvar1, mask, LAST, ierr )
         if (ierr /= 0) then
            call free_jac
            return
         end if

         call free_jac

   end subroutine semi_implicit_extrapolation

   subroutine free_jac
         implicit none
         call free_csr_matrix( reduced_jac_csr )
         call free_csc_matrix( reduced_jac_csc )
         call free_coo_matrix( reduced_jac_coo )
   end subroutine free_jac

   subroutine compute_error( k, y_k )
         integer, intent(in) :: k !< row of tableau
         integer :: idx
         real(r8), intent(in) :: y_k(nsp)
         real(r8) :: err_k(nsp)
         logical :: mask(nsp)

         mask     = y_k > 1e-27_r8
         err_k(:) = abs ( yerr(:,k) / yscal(:) )
         abuerr   = maxval( err_k(:), mask=mask ) / dgrd
         idx      = maxloc( err_k(:), dim = 1 )

         mask     = y_k < grdthreshold
         negerr   = maxval( -y_k, mask=mask ) / grdthreshold

         !> account for error in abundance sum, except for if the abundances
         !> have converged exceedingly well (to 1e-3 below the requested
         !> tolerance, e.g. 1e-6 when 1e-3 was requested)
         if (abuerr < 1e-3_r8) then
            !> abundance converged very well, so smaller time steps won't change
            !> much
            errmax = abuerr
         else
            !> take max of abundance error and mass conservation error
            sumerr   = abs(ONE - sortedsum(y_k)) / SUM_TOL
            errmax   = max(abuerr, sumerr)
         end if

         !> take max of negative abundance error and error so far
         errmax = max(errmax, negerr)

#ifdef VERBOSE
         !> reporting
         if (errmax == sumerr) then
            write(*,'(a5,3(a6,es15.7))') "bd-> ", " serr:", sumerr*SUM_TOL, " rerr:", sumerr, " yerr:", abuerr
         else if (errmax == negerr) then
            write(*,'(a5,3(a6,es15.7))') "bd-> ", " nerr:", negerr, " serr:", sumerr, " yerr:", abuerr
         else
            write(*,'(a8,a5,2x,4(es23.15))') 'bd err: ', zis(idx), yerr(idx,k), yscal(idx), abs( yerr(idx,k) / yscal(idx) / dgrd )
         end if
#endif
   end subroutine compute_error


   !> using Neville's algorithm, perform polynomial extrapolation of the resulting abundances from
   !> progressively higher order approximations to the solution
   subroutine extrapolate( k, yest, y_extrap )
         integer, intent(in) :: k !< row of tableau
         real(r8), intent(in) :: yest(nsp)
         real(r8), intent(out) :: y_extrap(nsp)
         integer :: &
               k1, & !< column of tableau
               j
         real(r8) :: delta, f1, f2, d(nsp), q, yz(nsp)

         x(k)           = xest
         yz(:)          = yest(:)
         yerr(:,k)      = yest(:)

         if ( k == 1 ) then
            ! first row, just a constant
            qcol(:,k)   = yest(:)
            y_extrap(:) = yest(:)
            return
         end if

         d(:) = yest(:)
         do k1 = 1, k - 1
            delta = ONE / ( x(k-k1) - xest ) 
            f1 = xest * delta
            f2 = x(k-k1) * delta
            do j = 1, nsp
               q          = qcol(j,k1)
               qcol(j,k1) = yerr(j,k)
               delta      = d(j)  - q
               yerr(j,k)  = f1    * delta
               d(j)       = f2    * delta
               yz(j)      = yz(j) + yerr(j,k)
            end do
         end do
         qcol(:,k)   = yerr(:,k)
         y_extrap(:) = yz(:)
   end subroutine extrapolate


   subroutine bader_deuflhard_step( y_k, delta_k, h, nvar1, mask, STAGE, ierr )
         use communication
         implicit none 
         real(r8) :: x(nsp), rhs(nsp), y_k(nsp), delta_k(nsp), h, dt_est
         real(r8), allocatable :: rhs_check(:), check(:)
         integer(i4) :: STAGE, nvar1, ierr, i, j, iopt, nrhs, ldb, info
         integer(kind=8) :: factors
         logical :: mask(nsp)

#ifdef USE_SUPERLU
         x(:) = ZERO; ierr = 0

         nrhs = 1; ldb = nvar1

         call calculate_dxdt( y_k, dt_est )

         select case ( STAGE )
         case( FIRST )
            rhs(:) = h * dxdt(:)
         case( MID, LAST )
            rhs(:) = h * dxdt(:) - delta_k(:)
         end select

         call reduce_rhs_bd( nvar1, mask, rhs )

         if ( STAGE == FIRST ) then
            ! factorize the matrix. The factors are stored in *factors* handle.
            iopt = 1
            call c_fortran_dgssv( iopt, nvar1, reduced_nnz, nrhs, &
               reduced_jac_csc%val, reduced_jac_csc%row, reduced_jac_csc%col, &
               reduced_rhs, ldb, factors, info )
         end if

         ! solve the system using the existing factors.
         iopt = 2
         call c_fortran_dgssv( iopt, nvar1, reduced_nnz, nrhs, &
            reduced_jac_csc%val, reduced_jac_csc%row, reduced_jac_csc%col, &
            reduced_rhs, ldb, factors, info )

         call ddecompac( mask, reduced_rhs, x, 1, nsp, nvar1, x )

         select case ( STAGE )
         case( FIRST, LAST )
            delta_k(:) = x(:)
         case( MID )
            delta_k(:) = delta_k(:) + TWO * x(:)
         end select

         y_k(:) = y_k(:) + delta_k

         if ( STAGE == LAST ) then
            iopt = 3
            call c_fortran_dgssv( iopt, nvar1, reduced_nnz, nrhs, &
                  reduced_jac_csc%val, reduced_jac_csc%row, reduced_jac_csc%col, &
                  reduced_rhs, ldb, factors, info )
         end if

         ! check for large abundances as an indication that the time step is too large
         if (any(y_k > TWO)) then
            ierr = HUGE_ABUNDS
            ! since the next time we'll be back here is with a different global time step,
            ! we should free the SuperLU storage
            if (STAGE /= LAST) then
               ! free the storage allocated inside SuperLU if we didnt already during LAST
               iopt = 3
               call c_fortran_dgssv( iopt, nvar1, reduced_nnz, nrhs, &
                     reduced_jac_csc%val, reduced_jac_csc%row, reduced_jac_csc%col, &
                     reduced_rhs, ldb, factors, info )
            end if
         end if

         ! free dynamic arrays
         call free_reduced_rhs()

#else
         stop "superLU is required for the bader-deuflhard solver"
#endif

   end subroutine bader_deuflhard_step


end module bader_deuflhard
