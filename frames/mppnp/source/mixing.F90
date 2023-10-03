module mixing

   ! for performing mixing in 1d models, either by:
   !     idifmix = 0: instantaneous mixing
   !     idifmix = 1

   use utils, only: r8, fmtprofile, wallclocktime
   use array_sizes, only: nsp, msl
   use frame_knobs, only: code_source, iolevel

   implicit none

   private

   ! choose instantaneous mixing approximation or diffusive mixing approximation
   integer, parameter :: idifmix = 1

   ! number of grid zones
   integer :: m
   ! lagrangian, eulerian diffusion coefficients, mass coordinate, abundances
   real(r8) :: d_lagr(msl), d_eul(msl), xm(msl), yps(msl,nsp)
   ! diffusion limiter (max. val. of lagrangian diffusion matrix element, right?)
   real(r8) :: sig_term_limit
   ! timing
   real(r8) :: tti1, tti2, mixing_time

   public do_mixing, sig_term_limit, mixing_time


contains



   subroutine do_mixing(xm, dq_reverse, rho, rad, d_eul_in, yps, m, dzeit)
         implicit none
         real(r8), dimension(msl) :: rho, rad, xm, dq_reverse, d_eul_in
         real(r8) :: yps(msl,nsp), dzeit
         integer :: m

         d_eul(:) = d_eul_in

         select case(idifmix)
         case(0)
            call mix_instantaneous(d_eul, m, xm, yps)
         case(1)
            tti1 = wallclocktime()
            call mix_calculate_D_Lagrangian(rho, rad, m)
            call mixdiff08(xm, dq_reverse, yps, m, dzeit)
            tti2 = wallclocktime()
            mixing_time = tti2 - tti1
            if (iolevel > 1) write(*,fmtprofile), 'mixdiff08 cpu_time/s = ',mixing_time
         case default
            stop 'unknown mixing choice for idifmix in mixing module'
         end select

   end subroutine do_mixing



   subroutine mix_calculate_D_Lagrangian(rho, rad, m)

         ! transform to Lagrangian coordinate: D_Lagr = (4*PI*R**2*RHO)**2*D_Eularian
         ! and also delta_xm in mixdiffnet is in solar mass units, and since there is 
         ! a second order spatial derivative:
         !     D_Lagr,mppn = (4*PI*R^2*RHO/Msun)**2*D_Eularian
         ! units from reading HPF files: [R]=Rsun, [rho]=cgs         
         !     4*Pi*Rsun^2/Msun = 3.0589703419947736e-11

         ! input:
         !     m:          number of mass shells
         !     d_eul:      eulerian diffusion coefficient
         !     rho:        density
         !     rad:        radius (in solar radii?)

         ! output:
         !     d_lagr:     Lagrangian diffusion coefficient

         integer :: j, m
         real(r8), parameter :: dconst = 3.0589703419947736e-11_r8
         real(r8) :: lnrho, lnrad
         real(r8), dimension(m), intent(in) :: rho, rad

         d_lagr(:) = 0._r8

         do j = 1, m
            if ( d_eul(j)  <=  0._r8 ) cycle
            if (j == 1) then
               lnrho = log(rho(j))
               lnrad = log(rad(j))
            else
               lnrho = (log(rho(j)) + log(rho(j-1))) / 2._r8
               lnrad = (log(rad(j)) + log(rad(j-1))) / 2._r8
            endif
            ! dppgLagr: Lagrangian coordinates; assuming solar mass units
            d_lagr(j) = d_eul(j) * dconst ** 2._r8 * exp(2._r8 * (2._r8 * lnrad + lnrho))
         end do

         return

   end subroutine mix_calculate_D_Lagrangian



   subroutine mixdiff08( xm, dq, yps, m, dzeit )

         ! $Id: mixdiffnet99.f,v 1.5 2001/04/02 21:20:30 fhg Exp fhg $
         !   Id: pde.f,v 1.2.1.1 1998/07/19 21:29:49 her Exp 
         ! *** solves the diffusion partial differential equation

         use frame_knobs
         use nuc_data

         implicit none

         integer :: m, j, insp

         real(r8) ::  a(msl), b(msl), c(msl) ! matrix coefficients: 
         real(r8) ::  as(msl), bs(msl), cs(msl) ! matrix coefficients: 
         ! b main diagonal
         real(r8) ::  yn(msl), ynp(msl),yps(msl,nsp)
         real(r8) ::  dzeit    ! stepsize delta_t in sec
         real(r8) ::  xm(msl), dq(msl) ! mass coordinate and shell mass

         !     *** Setting up Matrix: 
         if (code_source == 'GNV' .or. code_source == 'EVL') then
            call matrixdiffus(DZEIT, xm, d_lagr, a, b, c, m)
         else if(code_source == 'MES') then
            if ( .not. dq(1) > 0._r8 ) stop 'dq = 0. or NaN. Check your delta mass input.'
            call matrixdiffus_mesa(DZEIT, dq, d_lagr, a, b, c, m)
         end if

         diffuse: do insp=1,NSP-1
            if (considerisotope(insp)) then
               yn = yps(:,insp)
               as = a
               bs = b
               cs = c

               call tridag (as, bs, cs, yn, ynp, m)
               yps(:,insp) = ynp
            end if
         end do diffuse

         return

   end subroutine mixdiff08



   subroutine tridag( a, b, c, r, u, n )

         use array_sizes
         use utils, only: r8

         integer :: n, j
         real(r8) :: a(n), b(n), c(n), r(n), u(n)
         real(r8) :: bet, gam(MSL)

         if (b(1) == 0._r8) stop 'tridag: rewrite equations'
         bet = b(1)
         u(1) = r(1) / bet
         do j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            if (bet == 0._r8) stop 'Tridag failed.  The matrix is singular.'
            u(j) = (r(j)-a(j) * u(j-1)) / bet
         end do
         do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
         end do

         return

   end subroutine tridag



   subroutine  matrixdiffus( h, xm, D, a, b, c, m )

         ! compute coefficient matrix for diffusion solver
         ! method: fully implicit (see Numerical Recipes, p. 840)

         integer :: m, j
         real(r8) :: a(msl), b(msl), c(msl) ! matrix coefficients
         real(r8) :: d(msl)                 ! diffusion coefficent 
         real(r8) :: xm(msl)                ! star mass fraction
         real(r8) :: alpha,xmm,xr,xl       
         real(r8) :: h,xmdum(3) ! Mass coordinates for local gridpoints

         !     first row (center in PPN)    [siehe L0084, Eq.24,25]
         xmdum(2)=xm(1)
         xmdum(3)=xm(2)
         xr=xmdum(3)-xmdum(2)
         alpha=D(2)*h/(xr*xr)

         if (sig_term_limit  <  alpha) then
            print*, 'sig_term_limit reached by Dcoeff [D(2)*h/(xr*xr)]'
            print*, 'alpha:',2,alpha,'-->',sig_term_limit
            alpha = sig_term_limit
         endif

         a(1) = 0._r8
         b(1) = 1. + alpha
         c(1) = -1.d0*alpha

         ! all intermediate values
         do j=2,m-1
            xmdum(1)=xm(j-1)
            xmdum(2)=xm(j)
            xmdum(3)=xm(j+1)
            xl=xmdum(2)-xmdum(1)
            xr=xmdum(3)-xmdum(2)
            xmm=0.5D0 * (xmdum(3)-xmdum(1))
            alpha=h/xmm

            if (sig_term_limit  <  D(j+1)*alpha/xr) then
               if (iolevel  >=  2) then
                  print*, 'sig_term_limit reached by Dcoeff (D(j+1)*alpha/xr)'
                  print*, 'Dcoeff:',j+1,D(j+1),'-->',sig_term_limit/alpha*xr
                  print*, 'alpha:',alpha,', xr: ',xr,', sig: ',sig_term_limit
               endif
               D(j+1) = sig_term_limit/alpha*xr
            endif

            a(j)  =  -1.D0 * D(j) * alpha / xl
            b(j)  =   1.D0 + (D(j+1) / xr + D(j) / xl) * alpha
            c(j)  =  -1.D0 * D(j+1) * alpha / xr

            if (xl == 0.0_r8 .or. xr == 0.0_r8) then
               stop 'dm =0 in diffusion subroutine'
            endif

         end do

         ! last row (stellar surface in ppn)
         xmdum(1)=xm(M-1)
         xmdum(2)=xm(M)
         xl=xmdum(2)-xmdum(1)
         alpha=D(M)*h/(xl*xl)
         if (alpha >=  sig_term_limit) print*,'alphaM very large',alpha
         a(M) = -1._r8 * alpha
         b(M) = 1._r8 + alpha
         c(M) = 0._r8

         if (xl == 0._r8 .or. xr == 0._r8) then
            stop 'dm =0 in diffusion subroutine'
         endif

         return

   end subroutine matrixdiffus



   subroutine  matrixdiffus_mesa( h, dq, d, a, b, c, m )

         integer :: m, j
         real(r8) :: a(msl), b(msl), c(msl), & ! matrix coefficients
               d(msl), d_reduce(msl), & ! diffusion coefficient
               dq(msl), & ! delta mass
               xmm, xr, xl, & ! mass coordinate
               term_a,term_b,term_c,alpha,d_coeff, &
               h, xmdum(3) ! Mass coordinates for local gridpoints

         d_reduce(1)    = d(1)
         d_reduce(2:)   = 0._r8
         d_coeff        = 0._r8
         term_a         = 0._r8
         term_b         = 0._r8
         term_c         = 0._r8
         a              = 0._r8
         b              = 0._r8
         c              = 0._r8

         !     first row (center in PPN)    [siehe L0084, Eq.24,25]
         xr             = dq(1)
         alpha          = h / xr
         term_c         = 1._r8 / xr
         d_coeff        = d(2)

         call reduce_term(alpha, term_c, d_coeff)

         d_reduce(2)    = d_coeff

         a(1) =  0._r8
         b(1) =  1._r8 + term_c * alpha * d_reduce(2)         
         c(1) = -1._r8 * term_c * alpha * d_reduce(2)

         ! all intermediate values
         do j = 2, m - 1
            xl       = dq(j-1)
            xr       = dq(j)
            xmm      = 0.5_r8 * (dq(j) + dq(j-1))
            alpha    = h / xmm
            term_a   = 1._r8 / xl     ! D(j)/xl         
            term_c   = 1._r8 / xr     ! D(j+1)/xr
            d_coeff  = D(j+1)

            call reduce_term(alpha,term_c,d_coeff)

            d_reduce(j+1) = d_coeff

            term_b   = term_c * d_reduce(j+1) + term_a * d_reduce(j)

            a(j)     = -1._r8 * term_a * d_reduce(j) * alpha
            b(j)     =  1._r8 + term_b * alpha
            c(j)     = -1._r8 * term_c * d_reduce(j+1) * alpha

         end do

         !     last row (surface in PPN)    [siehe L0084, Eq.24,25]
         xl    = dq(m)
         alpha = h/(xl*xl)

         a(m)  = -1._r8 * alpha * d(m)
         b(m)  =  1._r8 + alpha * d(m)
         c(m)  =  0._r8

         return

   end subroutine matrixdiffus_mesa



   subroutine reduce_term(alpha, onedxr, d_coeff)

         real(r8) ::  alpha, onedxr, d_coeff, factor

         !*** dpa *** the parameter limit should not be too small, like 1.0d+2, otherwise dcoeff gets
         !            reduced to a level when chemical composition in both H and He convective cores
         !            is not uniform. With the value limit=1.0d+10, H and He convective cores are
         !            well mixed in one of test cases. We will work on this problem.

         if (alpha * onedxr * d_coeff > sig_term_limit) then
            d_coeff = sig_term_limit / (alpha * onedxr)
         end if

         return

   end subroutine reduce_term



   subroutine mix_instantaneous(d_eul, m, xm, yps)

         ! ^_^ FIXME: I think this has been broken for a while. the first two loops do not have i defined but it is used...

         !     kconvmax: max number of convective zones
         !     num_conv_zones: actual number of conv. zones
         !     kl(i): lower boundary (shell number) of conv. zone i
         !     ku(i): upper boundary (shell number) of conv. zone i
         !     xmconv: used to sum up mass in a conv. zone
         !     conv: used to calculate averaged abund. i na conv. zone

         ! input:
         !     m:       number of mass shells
         !     d_eul:   eulerian diffusion coefficient
         !     xm:      mass coordinate
         
         ! input/output:
         !     yps:     abundance (dimension msl,nsp) array

         integer, parameter :: kconvmax = 100
         integer :: i, j, m, num_conv_zones, kl(kconvmax), ku(kconvmax)
         real(r8) :: xmconv, yconv(nsp), d_eul(m), xm(msl), yps(msl,nsp)

         kl(:)    = -1
         ku(:)    = -1
         xmconv   = 0._r8
         yconv    = 0._r8
         num_conv_zones = 0
         
         ! ^_^ identify mesh points of convective boundaries
         if (d_eul(1) > 1.e10_r8) kl(1) = 1

         do j = 2, m
            ! lower boundary
            if (d_eul(j-1) < 1.e10_r8 .and. d_eul(j) > 1.e10_r8) then
               kl(i) = j
            endif
            ! upper boundary
            if (d_eul(j-1) > 1.e10_r8 .and. d_eul(j) < 1.e10_r8) then
               ku(i) = j
               ! count
               num_conv_zones = num_conv_zones + 1
            endif
         end do

         if ( kl(i)  >  0 .and. d_eul(m)  >  1.e10_r8 ) then
            ku(i) = m
            num_conv_zones = num_conv_zones + 1
         endif

         if (iolevel > 1) print *, 'num_conv_zones = ', num_conv_zones

         ! *** This loop goes over all defined convection zones and averages out
         !     the abundances
         do i = 1, num_conv_zones
            ! do not mix if conv. zone is only one cell wide
            if ( kl(i)  ==  ku(i) ) cycle
            yconv(:) = 0._r8
            xmconv = xm(ku(i)) - xm(kl(i))
            yconv(:) = yps(kl(i),:) * (xm(kl(i)+1) - xm(kl(i))) / 2._r8
            ! mass-weight
            do j = kl(i) + 1, ku(i) - 1
               yconv(:) = yconv(:) + yps(j,:) * (xm(j+1)-xm(j-1)) / 2._r8
            end do
            yconv(:) = yconv(:) + yps(ku(i),:) * (xm(ku(i)) - xm(ku(i)-1)) / 2._r8
            ! average
            yconv(:) = yconv(:) / xmconv
            do j = kl(i), ku(i)
               yps(j,:) = yconv(:)
            end do
         end do

         return

   end subroutine mix_instantaneous



end module mixing
