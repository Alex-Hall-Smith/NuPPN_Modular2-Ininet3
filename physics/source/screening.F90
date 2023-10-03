!>  Implement screening a la DeWitt, Graboske & Cooper 1973, ApJ, 181, 439 and
!>  Graboske, DeWitt, Grossman & Cooper 1973, ApJ, 1181, 457
!>
!>  @author Samuel Jones
!>
!>  Looking into this again, I see that we are assuming that the degeneracy
!>  is negligible, which is probably not the case in many of our situations.
!>  I have made this explicit by setting theta_e = 1. We should address this
!>  in the future. I also think that there is a factor of 1/abar missing in
!>  the definition of lambda0 (l0)

module screening
   use screening_graboske
   use screening_chugunov
   use physics_knobs
   use utils
   use constants

   implicit none
   private

   public screen_init, screen_precalc, calculate_screening

contains

   !> initialise screening module: calculate some constants and pre-process a
   !> hash table of proton numbers to the power of a few constants, that would
   !> be costly to calculate during the evolution loop of the code
   subroutine screen_init()
      implicit none

      select case(screen_option)
      case(NO_SCREENING)
         continue
      case(GRABOSKE)
         call screen_init_graboske()
      case(CHUGUNOV)
         call screen_init_chugunov()
      end select
   end subroutine screen_init


   subroutine screen_precalc(xin, tkin, rhoin)
         use nuc_data, only: considerisotope
      implicit none
      real(r8), intent(in) :: xin(nsp), tkin, rhoin

      select case(screen_option)
      case(NO_SCREENING)
         continue
      case(GRABOSKE)
         call screen_precalc_graboske(xin, tkin, rhoin)
      case(CHUGUNOV)
         call screen_precalc_chugunov(xin, tkin, rhoin, considerisotope)
      end select
   end subroutine screen_precalc


   subroutine calculate_screening( r1, r2, xn, tk, rho, screen )
      implicit none
      integer(i4) :: r1, r2 !< indices of reactants in network
      real(r8) :: z1, z2, a1, a2, x1, x2, xn(nsp), tk, rho, anetw, znetw, gam, screen
      common/cnetw/anetw(nsp),znetw(nsp)

      ! get charge numbers of reactants
      z1 = znetw(r1); z2 = znetw(r2)

      select case(screen_option)
      case(NO_SCREENING)
         screen = ONE
      case(GRABOSKE)
         call calculate_screening_graboske(z1, z2, screen)
      case(CHUGUNOV)
         call calculate_screening_chugunov(z1, z2, gam, screen)
      end select
      
      ! it is with a heavy heart that I still need to limit the screening function
      screen = min(screen, 1.e30_r8)
   end subroutine calculate_screening

end module screening
