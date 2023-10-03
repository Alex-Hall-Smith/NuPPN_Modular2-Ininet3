!> @Author Samuel Jones

! verbose mode
#undef VERBOSE

module nse_solver
   use array_sizes
   use nuc_data
   use rates
   use utils
   use constants
   use solver_knobs
   use physics_knobs
   use frame_knobs
   use solver_diagnostics
   use errors
   use backward_euler, only: newton_iteration
   use nse_torch, only: testnse
   use nse_swj, only: compute_nse, nse_num_isos
   use communication
   implicit none
   private

   real(r8) :: lm_ov_a(nsp), anetw(nsp), znetw(nsp)

   common /cnetw/ anetw, znetw                                

   public nse_solver_cache_lm_ov_a, nse_yedot_rk4, nse_yedot_euler

contains

   !> calculate and save \f$ \frac{1}{A}(\lambda_i^\beta - \lambda_i^\mathrm{ec}) \f$ for each nuclide
   subroutine nse_solver_cache_lm_ov_a()
         use nuc_data, only: considerreaction
         use rates, only: v
        use reaction_info, only: ilabb, i_bm, i_ec, i_ve, i_vp, i_vsp, i_vsn, i_vsa
         integer(i4) :: i, par, dau
         real(r8) :: lm, zp, zd, a

         ! set all the rates that are not weak rates to zero
         do i = 1, nre
            select case(ilabb(i))
            case(i_bm,i_ec,i_ve,i_vp)
               ! weak reaction - keep this one
               continue
            case default
               v(i) = ZERO
            end select
         end do

         lm_ov_a(:) = ZERO

         do i = 1, nre
            if ( .not. considerreaction(i) ) cycle
            if ( v(i) == ZERO ) cycle
            lm = v(i); par = k1(i); dau = k7(i)
            zp = znetw(par); zd = znetw(dau); a = anetw(par)

            if ( zp < zd ) then
               ! beta-decay or positron-capture
               lm_ov_a(par) = lm_ov_a(par) + lm / a
            else
               ! electron-capture or positron-decay
               lm_ov_a(par) = lm_ov_a(par) - lm / a
            end if
         end do
   end subroutine nse_solver_cache_lm_ov_a

   !> Integrate nuclear network in time assuming nuclear statistical equilibrium and
   !> coupling non-equilibrium weak rates using RK4 time integrationi of the equation
   !> \f$ \dot{Y}_\mathrm{e} = \sum_i X_i(Y_\mathrm{e}) \cdot (\lambda_i^\beta - \lambda_i^\mathrm{ec})
   !>
   !> This routine performs a Cash-Karp 6-step 5th order accurate Runge Kutta integration, which
   !> gives a 4th-order solution with error estimate, but we use the 5th order solution
   subroutine nse_yedot_rk4( nvar, yps, t9, rho, ye, dt, dtrcmd, ierr )
         real(r8), dimension(nsp) :: yps, y1nse
         real(r8) :: dt, ye, ye_old, rho, t9, mu_p, mu_n, error, fac, dtrcmd
         real(r8), parameter :: &
               maxfac = 2._r8, &  !< maximum factor to recommend increasing time step by
               minfac = 0.1_r8    !< minimum factor to recommend reducing time step by
         integer(i4) :: nvar, ierr, iterations

         ye_old = ye ; y1nse(:) = yps(:); ierr = 0

         ! RK5 step and error (for 4th order rk):
         call rk4_nse_step( dt, rho, t9, ye, y1nse, nvar, error, ierr )
         if ( ierr /= 0 ) then
            error = ye_old
         end if

         ! check error
         if ( ierr /= 0 .or. ye < ZERO .or. ye > ONE ) then
            ! singular NSE jacobian matrix or time step too large that unphysical Y_e attained
            if ( ye < ZERO .or. ye > ONE ) then
               ierr = BAD_YE
            end if
            ierr = BAD_YE
            fac = HALF * abs( ye_tol * ye_old / max(error,SMALL) ) ** FIFTH
         else if ( error / ye > ye_tol ) then
            ! poor convergence
            ierr = NSE_WEAK
            fac = HALF * abs( ye_tol * ye / max(error,SMALL) ) ** FIFTH
         else
            ! solution is good
            ierr = 0
            fac = abs( ye_tol * ye / max(error,SMALL) ) ** FIFTH
            yps(:) = y1nse(:)
         end if
         dtrcmd = dt * max(min(fac,maxfac),minfac)
   end subroutine nse_yedot_rk4


   function yedot( x )
         use nuc_data, only: considerisotope
         real(r8) :: x(nsp), yedot
         integer :: i
         yedot = sortedsum( x(:) * lm_ov_a(:), mask = considerisotope )

#ifdef VERBOSE
         print *, 'yedot = ', yedot
#endif
   end function


   subroutine rk4_nse_step( dt, rho, t9, ye, y1nse, nvar, error, ierr )
         real(r8), dimension(nsp) :: y1nse
         real(r8) :: dt, ye, ye_0, rho, t9, mu_p, mu_n, &
               rk_k1, rk_k2, rk_k3, rk_k4, rk_k5, rk_k6, ye_5th, ye_4th, error, tti1
         !> RK4/5 coefficients
         real(r8), parameter :: &
               b21 = FIFTH, &
               b31 = 3._r8 / 40._r8, &
               b32 = 9._r8 / 40._r8, &
               b41 = 3._r8 / 10._r8, &
               b42 = -9._r8 / 10._r8, &
               b43 = SIX_FIFTHS, &
               b51 = -11._r8 / 54._r8, &
               b52 = 5._r8 * HALF, &
               b53 = -70._r8 / 27._r8, &
               b54 = 35._r8 / 27._r8, &
               b61 = 1631._r8 / 55296._r8, &
               b62 = 175_r8 / 512._r8, &
               b63 = 575._r8 / 13824._r8, &
               b64 = 44275._r8 / 110592._r8, &
               b65 = 253._r8 / 4096._r8, &
               c1  = 37._r8 / 378._r8,   c1s = 2825._r8 / 27648._r8, &
               c2  = ZERO,               c2s = ZERO, &
               c3  = 250._r8 / 621._r8,  c3s = 18575._r8 / 48384._r8, &
               c4  = 125._r8 / 594._r8,  c4s = 13525._r8 / 55296._r8, &
               c5  = ZERO,               c5s = 277._r8 / 14336._r8, &
               c6  = 512._r8 / 1771._r8, c6s = FOURTH
         integer :: ierr, iterations, nvar

         ye_0 = ye
         ! ^_^ calculate y1neu NSE distribution at this t9, rho, ye
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return

         ! ^_^ RK4
         rk_k1 = dt * yedot( y1nse )

         ye = ye_0 + ( b21 * rk_k1 )
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return
         rk_k2 = dt * yedot( y1nse )

         ye = ye_0 + ( b31 * rk_k1 + b32 * rk_k2 )
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return
         rk_k3 = dt * yedot( y1nse )

         ye = ye_0 + ( b41 * rk_k1 + b42 * rk_k2 + b43 * rk_k3 )
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return
         rk_k4 = dt * yedot( y1nse )

         ye = ye_0 + ( b51 * rk_k1 + b52 * rk_k2 + b53 * rk_k3 + b54 * rk_k4 )
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return
         rk_k5 = dt * yedot( y1nse )

         ye = ye_0 + ( b61 * rk_k1 + b62 * rk_k2 + b63 * rk_k3 + b64 * rk_k4 + b65 * rk_k5 )
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return
         rk_k6 = dt * yedot( y1nse )

         ! 5th order and 4th order solutions
         ye_5th = ye_0 + c1  * rk_k1 + c2  * rk_k2 + c3  * rk_k3 + c4  * rk_k4 + c5  * rk_k5 + c6  * rk_k6
         ye_4th = ye_0 + c1s * rk_k1 + c2s * rk_k2 + c3s * rk_k3 + c4s * rk_k4 + c5s * rk_k5 + c6s * rk_k6

         ! get new abundances at 5th-order ye solution:
         ye = ye_5th
         select case ( nse_option )
         case ( NSE_TORCH )
            call testnse(t9, rho, ye, nvar, y1nse, ierr)
         case ( NSE_SWJ )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         if ( ierr /= 0 ) return

         ! error estimate (actually an estimate of error in 4th-order solution)
         error = abs( ye_5th - ye_4th )
#ifdef VERBOSE
         print *, 'NSE: RK4 weak couplling error: ', error * ye
#endif
   end subroutine rk4_nse_step


   subroutine nse_yedot_euler( nvar, yps, t9, rho, ye, dt, dtrcmd, ierr )
         real(r8), dimension(nsp) :: yps, y1nse
         real(r8) :: dzeit, deltamax, ye, d_ye, dyedt, lm, zp, &
               zd, dt, tdone, a, term, rho, t9, mu_p, mu_n, dtrcmd
         integer :: nvar, i, ierr, iterations
         logical :: done

         ! allowing the Ye to only change by 10^-5 per subtimestep shows
         ! converged results
         deltamax = 1.e-5_r8

         ! calculate y1neu NSE distribution at this t9, rho, ye
         select case ( nse_option )
         case( 0 )
            call testnse(T9, rho, ye, nvar, y1nse, ierr)
         case( 1 )
            call compute_nse(t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr)
         end select
         ! check error
         if ( ierr /= 0 .or. ye < ZERO .or. ye > ONE ) then
            ! singular NSE jacobian matrix or time step too large that unphysical Y_e attained
            ierr = BAD_STEP
            dtrcmd = 0.1_r8 * dt
            return
         end if

         dyedt = yedot( y1nse )
         d_ye  = dyedt * dt

         if ( d_ye > deltamax ) then
            ierr = NSE_WEAK
            dtrcmd = deltamax / abs( dyedt )
            return
         else
            ye = ye + d_ye
         end if

         select case( nse_option )
         case( NSE_TORCH )
            call testnse(T9, rho, ye, nvar, y1nse, ierr)
         case( NSE_SWJ )
            call compute_nse( t9, rho, ye, y1nse, mu_p, mu_n, iterations, ierr )
         end select
         
         select case( ierr )
         case ( 0 )
            yps(:) = y1nse(:)
         case default
            ierr = NSE_CNVRG
         end select
         dtrcmd = deltamax / abs( dyedt )

   end subroutine nse_yedot_euler

end module nse_solver
