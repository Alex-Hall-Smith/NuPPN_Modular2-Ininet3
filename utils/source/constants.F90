!> @author S. W. Jones
!
!> @brief Numerical and physical constants
!> @details This module provides numerical and physical constants to any other module or routine in the code that should need them

module constants
   use utils, only: r8, i4
   implicit none

   ! conversions
   real(r8), parameter :: &
         ev2erg   = 1.602177e-12_r8, &       !< conversion eV to erg
         kev2erg  = ev2erg * 1.0e3_r8, &     !< conversion keV to erg
         mev2erg  = ev2erg * 1.0e6_r8, &     !< conversion MeV to erg
         erg2mev  = 1._r8 / mev2erg, &       !< conversion erg to MeV
         mus2sec  = 1.e-6_r8, &              !< conversion \f$ \mu s \f$ to seconds
         msec2sec = 1.e-3_r8, &              !< conversion ms to seconds
         min2sec  = 60._r8, &                !< conversion minutes to seconds
         hrs2sec  = min2sec * 60._r8, &      !< conversion hours to seconds
         day2sec  = hrs2sec * 24._r8, &      !< conversion days to seconds
         yrs2sec  = day2sec * 365._r8        !< conversion years to seconds

   ! physical constants (in CGS units)
   real(r8), parameter :: &
         pi          = acos(-1._r8), &
         twopi       = 2._r8 * pi, &
         clight      = 2.99792458e10_r8, &             !< speed of light in vacuum (cm/s)
         clight2     = clight * clight, &
         avogadro    = 6.022140857e23_r8, &            !< Avogadro's number \f$ (\mathrm{mol}^{-1}) \f$
         boltzmann   = 1.3806485279e-16_r8, &          !< Boltzmann's constant \f$ (\mathrm{erg K}^{-1}) \f$
         deltap      = 7.288969_r8, &                  !< proton mass excess (MeV) from Audi & Wapstra 1995
         deltan      = 8.071323_r8, &                  !< neutron mass excess (MeV) from Audi & Wapstra 1995
         planck_mev  = 4.135667662e-21_r8, &           !< planck's constant (MeV s)
         planck_erg  = planck_mev * mev2erg, &         !< planck's constant (erg s)
         hbar_mev    = planck_mev / 2._r8 / pi, &      !< planck's constant over 2pi (MeV s)
         hbar_erg    = planck_erg / 2._r8 / pi, &      !< planck's constant over 2pi (MeV s)
         amu         = 1.6605402e-24_r8, &             !< atomic mass unit (g)
         mp_mev      = 938.27231_r8, &                 !< proton rest mass energy equivalence (MeV)
         mn_mev      = 939.56563_r8, &                 !< neutron rest mass energy equivalence (MeV)
         me_mev      = 0.5109989461_r8, &              !< electron rest mass energy equivalence (MeV)
         mp_g        = mp_mev * mev2erg / clight2, &   !< proton rest mass (g)
         mn_g        = mn_mev * mev2erg / clight2, &   !< neutron rest mass (g)
         me_g        = 9.10938356e-28_r8, &            !< electron rest mass (g)
         qe          = 4.80320425e-10_r8, &            !< electron charge (statcoulombs)
         qe2         = qe*qe                           !< electron charge**2 (statcoulombs**2)

   ! computer numbers
   integer(i4), parameter :: &
         r8size = 8, &
         i4size = 4, &
         ch5size = 5


   ! numbers and fractions
   real(r8), parameter :: &
         ZERO           = 0._r8, &
         ONE            = 1._r8 ,&
         TWO            = 2._r8, &
         THREE          = 3._r8, &
         FOUR           = 4._r8, &
         TEN            = 10._r8, &
         HALF           = 0.5_r8, &
         THIRD          = 1._r8 / 3._r8, &
         FOURTH         = 1._r8 / 4._r8, &
         FIFTH          = 1._r8 / 5._r8, &
         SIXTH          = 1._r8 / 6._r8, &
         SEVENTH        = 1._r8 / 7._r8, &
         EIGTH          = 1._r8 / 8._r8, &
         TWO_THIRDS     = 2._r8 * THIRD, &
         FOUR_THIRDS    = 4._r8 * THIRD, &
         FIVE_THIRDS    = 5._r8 * THIRD, &
         TWO_FIFTHS     = 2._r8 * FIFTH, &
         THREE_FIFTHS   = 3._r8 * FIFTH, &
         FOUR_FIFTHS    = 4._r8 * FIFTH, &
         SIX_FIFTHS     = 6._r8 * FIFTH, &
         FIVE_SIXTHS    = 5._r8 * SIXTH, &
         TWO_SEVENTHS   = 2._r8 * SEVENTH, &
         THREE_EIGTHS   = 3._r8 * EIGTH, &
         ln2            = log(TWO), &
         ln10           = log(10._r8), &
         SMALL          = 1.e-30_r8

end module constants
