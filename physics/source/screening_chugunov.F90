!>  Implement screening a la Chugunov, DeWitt & Yakovlev 2007, PhRvD, 76, 025028
!>
!>  @author Samuel Jones
!>

module screening_chugunov
   use utils
   use constants

   implicit none
   private

   real(r8) :: &
         fourov3pi, &      !< constant 4/(3*pi)
         fourpiqe2, &      !< constant 4*pi*qe*qe
         bohrconst, &      !< constant part of the bohr radius formula
         zetalim, &        !< upper lim of zeta corresponding to min. T/Tp ratio
         allzpow13(200), & !< the cube root of numbers 1-200
         abar, &           !< mean atomic weight
         zbar, &              !< mean charge
         a_av, &       !< average ion sphere radii
         a_e, &        !< electron sphere radius
         ntot, &       !< total ion number density
         tp, wp, &     !< plasma temperature and frequency
         tk, &         !< temperature to use for screening (may be corrected for limits)
         tn, &         !< normalised temperature (tk/tp)
         gamfac, &     !< factor of gamma, the coulomb coupling parameter, for pre-processing
         zeta          !< function of plasma temperature

   ! there are various realms of applicability in the Chugunov paper, for
   ! the particle-in-cell Monte-Carlo calculations, for the MKB
   ! approximation, and for the Potekhin-Cabrier-based fit to the MKB
   ! results for the mean field potential H(r).
   ! hence, here are parameters for smoothly limiting the domain of application
   real(r8), parameter :: &
      xlim = 1.e-99_r8, &    !< lower limit for mass fraction of reactant for screening
      gamfitlim = 600._r8, & !< upper limit for applicability of fit for h(gamma)
      g0 = 599._r8, &        !< centre value of sigmoid function for gammatilda
      t2 = 0.2_r8, &         !< minimum allowed ratio T/T_p before fading out
      t1 = 0.1_r8, &         !< floor value of T/T_p (T_p = plasma temperature)
      t0 = (t1 + t2)*HALF    !< centre value of sigmoid fade-out function for temperature

   logical :: &
         screening_initialised = .false.

   public screen_init_chugunov, calculate_screening_chugunov, screen_precalc_chugunov

contains

   !> initialise screening module: calculate some constants and pre-process a
   !> hash table of proton numbers to the power of one third
   subroutine screen_init_chugunov()
         implicit none
         integer(i4) :: i

         fourov3pi = FOUR/(THREE*pi*pi)
         fourpiqe2 = FOUR*pi*qe2
         bohrconst = hbar_erg**2/qe2

         do i = 1, 200
            allzpow13(i) = real(i,r8) ** THIRD
         end do

         screening_initialised = .true.
   end subroutine screen_init_chugunov

   subroutine screen_precalc_chugunov(x_in, tkin, rho, mask)
         use utils, only: calculate_azbar
         implicit none
         real(r8) :: &
               x_in(nsp), &  !< mass fractions
               anetw(nsp), & !< aion
               znetw(nsp), & !< zion
               rho, &        !< mass density (g/cc)
               tkin, &       !< temperature (K)
               s, &          !< evaluated sigmoid function
               denom
         logical :: mask(nsp)  !< which isotopes to consider
         common/cnetw/anetw,znetw

         !> calculate abar and zbar
         call calculate_azbar(x_in, anetw, znetw, abar, zbar, mask)

         !> calculate effective ion and electron sphere radii using abar and zbar of
         !> the multicomponent plasma, see Itoh et al. 1979 eqs. (1-3)
         ntot  = rho / (abar*amu) !< total ion number density
         denom = (THREE/(FOUR*pi*zbar*ntot))**THIRD
         a_av  = denom * zbar**THIRD
         a_e   = denom

         ! plasma frequency and temperature (chugunov 2007 eq 2)

         wp = (fourpiqe2*zbar*zbar*ntot/(abar*amu))**HALF
         tp = hbar_erg*wp/boltzmann
         tn = tkin/tp

         ! smooth cutoff for normalised temperature (T/T_p) between two points
         ! because this is an exponential call, try to limit the number of times it
         ! is used by only calling when we're sure tn is such that the sigmoid
         ! function is probably not equal to 1. To determine that we just use twice the
         ! upper limit of the fade region (2*t2) and half the lower limit (.5*t1)

         if (tn >= TWO*t2) then

            ! well above the blend, so sigmoid = 1

            tk = tkin

         else if (tn <= HALF*t1) then

            ! well below the blend, so sigmoid = 0

            tk = t1*tp

         else

            ! perform the fade 

            s = sigt(tn)
            tn = (ONE - s)*t1 + s*tn

            ! define the absolute temperature based on the faded normalised temperature

            tk = tn*tp

         end if

         gamfac = qe2/(a_av*tk)/boltzmann

         ! zeta (chugunov 2007 eq 3)

         zeta = (fourov3pi*tp*tp/tk/tk)**THIRD

   end subroutine screen_precalc_chugunov

   !> TODO: According to the end of Itoh et al. 1979, to generalise to a multi-component
   !> plasma the ion and electron sphere radii should be computed using Z = zbar, A = zbar,
   !> etc, but the gammas using z1 and z2 (the charges of the reactants)
   !> the calculation of abar and zbar should be done in a precalc routine
   subroutine calculate_screening_chugunov(z1, z2, gam, screen)
      implicit none
      ! coefficients from chugunov
      real(r8), parameter :: &
         c_a1 = 2.7822_r8, &
         c_a2 = 98.34_r8, &
         c_a3 = sqrt(THREE) - c_a1/sqrt(c_a2), &
         c_b1 = -1.7476_r8, &
         c_b2 = 66.07_r8, &
         c_b3 = 1.12_r8, &
         c_b4 = 65._r8, &
         alfa = 0.022_r8
      real(r8) :: &
         z1, z2, &     !< charge numbers of reactants
         tkin, &       !< temperature (K) input
         gam, &        !< on return, coulomb coupling parameter gamma
         screen, &     !< on return, screening factor for this reaction
         z1z2, &       !< z1*z2
         t1, t2, t3, & !< the three terms in the fitting formula for h0fit
         h0fit, &      !< mean field fit
         denom, &
         beta, gama, & !< coeffs of zeta formula
         gamtild, &    !< "effective" coulomb coupling parameter
         gamtild2, &   !< "effective" coulomb coupling parameter squared
         s             !< evaluated sigmoid function

      ! check for initialisation
      if (.not. screening_initialised) then
         stop "screening_chugunov: not initialised"
      end if

      ! check whether both reactants are charged ions

      if (z1 <= ZERO .or. z2 <= ZERO) then
         screen = ONE
         return
      end if

      z1z2 = z1*z2

      ! coulomb coupling parameter gamma (itoh 1979 eq 4)
      
      gam = z1z2*gamfac

      ! coefficients of zeta dependent on the ion coulomb coupling
      ! parameter (gam) (chugunov 2007, just after eq 21)

      beta = 0.41_r8 - 0.6_r8/gam
      gama = 0.06_r8 + 2.2_r8/gam

      ! gamma tilda (chugunov 2007 eq 21)

      gamtild  = gam/(ONE + alfa*zeta + beta*zeta**2 + gama*zeta**3)**THIRD

      ! smoothly cap the gammatilda value to 600 using sigmoid.
      ! again only do if we really need to, to avoid unnecessary exponential
      ! function evaluations

      if (gamtild >= 650) then

         ! sigmoid should be 1

         gamtild = gamfitlim

      else if (gamtild <= 550) then

         ! sigmoid should be 0, so keep gammatild as is

         continue
      else

         ! smoothly limit gammatilda using sigmoid function blend

         s = sigg(gamtild)
         gamtild = s*gamfitlim + (1-s)*gamtild

      end if

      gamtild2 = gamtild*gamtild

      ! for for mean field potential H (chugunov 2007 eq 19)

      t1 = gamtild**(THREE*HALF)*( c_a1/sqrt(c_a2+gamtild) + c_a3/(1 + gamtild) )
      t2 = c_b1*gamtild2 / (c_b2 + gamtild)
      t3 = c_b3*gamtild2 / (c_b4 + gamtild2)
      h0fit = t1 + t2 + t3
      screen = exp(h0fit)
      
      ! finally, limit screening factor to 1e30, warning if this was necessary
      screen = min(screen, 1e30_r8)
      !if (screen == 1e300_r8) then
         !print *, "warning: screening_chugunov: large screening factor was clipped to 1e300"
      !end if

   end subroutine calculate_screening_chugunov

   !> sigmoid function for fading out temperature
   real(r8) function sigt(x)
      implicit none
      real(r8) :: x

      sigt = ONE / (ONE + exp(-(x - t0)/0.01_r8))
   end function sigt

   !> sigmoid function for fading out gammatilda
   real(r8) function sigg(g)
      implicit none
      real(r8) :: g

      sigg = ONE / (ONE + exp(-(g - g0)/ONE))
   end function sigg

end module screening_chugunov



