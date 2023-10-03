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

module screening_graboske
   use array_sizes, only: nsp
   use nuc_data
   use utils
   use constants

   implicit none
   private

   real(r8) :: x(nsp)   ! mass fractions
   real(r8) :: tk       ! temperature (K)
   real(r8) :: rho      ! density (g/cc)

   ! various constants, parameters and products thereof for the screening
   ! calculations
   real(r8) :: f22
   real(r8) :: ztilda, ztildal0, l0, zbar, z2bar, zzmui
   real(r8) :: &
         const01, const02, const03
   real(r8) :: h120, h120is, h120is_pre
   real(r8) :: theta_e = 1.0_r8
   real(r8), parameter :: tk_screen_limit = 1.e6_r8
   real(r8) :: allz(200), allzpow23(200), allzpow43(200), &
         allzpow53(200), allzpow1p86(200)

   logical :: screening_initialised = .false.

   real(r8) :: an(nsp), zn(nsp)
   common/cnetw/an, zn

   public screen_init_graboske, screen_precalc_graboske, &
      calculate_screening_graboske

contains

   subroutine screen_init_graboske()
         ! ^_^ initialise screening module: calculate some constants and pre-process a
         !     hash table of proton numbers to the power of a few constants, that would
         !     be costly to calculate during the evolution loop of the code

         integer :: i

         f22         = 0.624_r8

         ! ^_^ pre-process exponents of z
         do i = 1, 200
            allz(i) = dble(i)
         end do
         allzpow23(:)   = allz(:) ** TWO_THIRDS
         allzpow43(:)   = allz(:) ** FOUR_THIRDS
         allzpow53(:)   = allz(:) ** FIVE_THIRDS
         allzpow1p86(:) = allz(:) ** 1.86_r8

         screening_initialised = .true.
   end subroutine screen_init_graboske



   subroutine screen_precalc_graboske(xin, tkin, rhoin)
         implicit none
         real(r8), intent(in) :: xin(nsp), tkin, rhoin
         real(r8)             :: xdiva
         integer              :: i

         ! check for initialisation

         if (.not. screening_initialised) then
            stop "screening_graboske: not initialized"
         end if

         ! init
         x(:)   = xin(:); tk       = tkin; rho   = rhoin
         zbar   = ZERO; z2bar      = ZERO; zzmui = ZERO; ztilda   = ZERO; h120 = ZERO;
         h120is = ZERO; h120is_pre = ZERO; l0    = ZERO; ztildal0 = ZERO

         if (minval(x) < -1.e-18_r8) then
            stop 'in screening: large negative abund. something is wrong. god help you'
         end if

         ! for strong screening
         do i = 1, nsp
            if ( an(i) > ZERO ) then
               xdiva      = x(i) / an(i)

               ! DeWitt, Graboske & Cooper 1973 eq 4b,c
               zbar       = zbar    + ( zn(i) * xdiva )

               zzmui      = zzmui   + xdiva
               z2bar      = z2bar   + ( zn(i) * zn(i) * xdiva )
               h120is_pre = h120is_pre + ( zn(i) ** 1.58_r8 * xdiva ) ! (for weak)
            end if
         end do

         ! DeWitt, Graboske & Cooper 1973 eq 4a
         ztilda      = (z2bar + zbar * theta_e) ** HALF

         ! DeWitt, Graboske & Cooper 1973 eq 6
         l0          = 1.88e8_r8 * ( rho / ( tk * tk * tk ) ) ** HALF
         ztildal0    = ztilda * l0

         ! Graboske, DeWitt, Grossman & Cooper 1973, Table 4
         const01     = 0.316_r8 * (zbar / zzmui) ** THIRD
         const02     = 0.737_r8 / (zbar * (l0 / zzmui) ** TWO_THIRDS)
         const03     = f22 * zbar ** THIRD * l0 ** TWO_THIRDS

         ! for weak screening
         h120is_pre = h120is_pre / (ztilda ** 0.58_r8 * zbar ** 0.28_r8)
         h120is_pre = h120is_pre * 0.38_r8 * l0 ** 0.86_r8
   end subroutine screen_precalc_graboske


   subroutine calculate_screening_graboske( z1, z2, screen )
         real(r8) :: z1, z2, z1pz2, screen
         real(r8) :: l12, h120, h120a
         !  pre-processed exponents:
         integer  :: iz1, iz2, iz1pz2
         real(r8) :: z1pow23, z1pow43, z1pow53, z2pow23, z2pow43, z2pow53, &
               z1pz2pow23, z1pz2pow43, z1pz2pow53, z1pow1p86, z2pow1p86, z1pz2pow1p86

         screen = ONE

         if ( tk < tk_screen_limit ) then
            return
         end if

         if (z1 <= ZERO .or. z2 <= ZERO) then
            return
         end if

         z1pz2 = z1 + z2

         ! DeWitt, Graboske & Cooper 1973 eq 6
         l12   = z1 * z2 * ztildal0

         iz1 = int(z1); iz2 = int(z2); iz1pz2 = int(z1pz2)
         h120 = ZERO

         if (l12 < 0.1_r8) then
            ! ^_^ apply weak screening
            h120        = l12
            screen      = exp(h120)
         else
            ! ^_^ calculate intermediate screening  using pre-processed exponents
            z1pz2pow1p86   = allzpow1p86(iz1pz2)
            z1pow1p86      = allzpow1p86(iz1)
            z2pow1p86      = allzpow1p86(iz2)
            h120is         = h120is_pre * (z1pz2pow1p86 - z1pow1p86 - z2pow1p86)
            if ( l12 <= 2._r8 ) then
               ! ^_^ apply intermediate screening
               screen      = exp(h120is)
            else
               ! ^_^ calculate strong screening without pre-processed exponents
               !            h120a =           ( z1pz2**(5/3) - z1**(5/3) - z2**(5/3) ) + &
               !                    const01 * ( z1pz2**(4/3) - z1**(4/3) - z2**(4/3) ) + &
               !                    const02 * ( z1pz2**(2/3) - z1**(2/3) - z2**(2/3) )

               ! ^_^ calculate strong screening using pre-processed exponents
               z1pow23 = allzpow23(iz1); z2pow23 = allzpow23(iz2)
               z1pow43 = allzpow43(iz1); z2pow43 = allzpow43(iz2)
               z1pow53 = allzpow53(iz1); z2pow53 = allzpow53(iz2)

               z1pz2pow23 = allzpow23(iz1pz2)
               z1pz2pow43 = allzpow43(iz1pz2)
               z1pz2pow53 = allzpow53(iz1pz2)

               h120a = ( z1pz2pow53 - z1pow53 - z2pow53 ) + &
                     const01 * ( z1pz2pow43 - z1pow43 - z2pow43 ) + &
                     const02 * ( z1pz2pow23 - z1pow23 - z2pow23 )
               h120a    = h120a * const03

               ! ^_^ apply either min or only strong, depending on value of Lambda12
               if (l12 <= 5._r8) then
                  if (h120is <= h120a) then
                     h120     = h120is
                  else
                     h120     = h120a
                  end if
               else if (l12 > 5._r8) then
                  h120 = h120a
               end if

               screen      = exp(h120)
            end if
         end if

   end subroutine calculate_screening_graboske


end module screening_graboske
