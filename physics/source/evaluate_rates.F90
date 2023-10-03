#undef VERBOSE
#undef PROFILE
module evaluate_rates
   use utils, only: r8, i4, wallclocktime, fmtprofile
   use array_sizes, only: nsp, nre
   use communication
   implicit none
   private

   public :: evaluate_all_rates

contains

   subroutine evaluate_all_rates( ye, nvar, nvrel, rho, t9, yps, nu )
      use frame_knobs, only: iolevel
      use physics, only: rnetw2008
      use reverse, only: calculate_reverse_rates
      use physics_knobs, only: detailed_balance, weak_rates, strong_rates, decay
      use vital, only: vital_rates_derivs
      use screening, only: screen_precalc, calculate_screening
      use rates, only: ant, anf, znt, znf, v, k1, k2, k3
      use nuc_data, only: ispe, considerreaction, zis
      use reaction_info
      use neutrinos
      implicit none
      real(r8), intent(in) :: ye, rho, t9, yps(nsp)
      integer(i4), intent(in) :: nvar, nvrel
      integer(i4) :: &
         r1, r2, & ! reactants
         i, &
         ihe4, ibe8 !< indices of he4 and be8 in the network, for special cases
      real(r8) :: &
         an(nsp), zn(nsp), Fscreen, Fscreenaux, tk, &
         tti1, tti2, rnetw2008_time, vital_time, reverse_time
      type(neutrino) :: nu
      common / cnetw / an, zn

      tti1 = wallclocktime()
      call rnetw2008(ye, an, zn, nvar, nvrel, rho, t9, yps, nu)
      rnetw2008_time = wallclocktime() - tti1

#ifdef VERBOSE
      write( *, fmtprofile ) 'rnetw2008 cpu_time/s = ', rnetw2008_time
#endif

      call apply_rate_factors()

      if (detailed_balance) then
         tti1 = wallclocktime()
         call calculate_reverse_rates(rho, t9, ye)
         reverse_time = wallclocktime() - tti1
#ifdef VERBOSE
         write( *, fmtprofile ) 'reverse   cpu_time/s = ', reverse_time
#endif
      end if

      tti1 = wallclocktime()

      ! calculate and apply screening factors
      tk = t9*1.e9_r8
      call screen_precalc(yps, tk, rho)

      ! locate helium 4 and berillium 8 for special case reactions
      ihe4 = ispe("HE  4"); ibe8 = ispe("BE  8")

      !> loop over all the reactions and screen them
      !> note that in previous versions, the VITAL routine made its own
      !> screening calls
      do i = 1, nre

         ! screen VITAL reactions whatever the conditions

         if (.not. considerreaction(i) .and. i > 110) cycle

         ! don't screen reactions with teeny rates

         if (v(i) <= 1.e-99_r8) cycle

         ! special cases in vital: 
         ! * tripla alpha
         ! * n(2a,g)be9

         select case(i)
         case(43)

            ! triple alpha

            call calculate_screening(ihe4, ihe4, yps, tk, rho,   Fscreen)
            call calculate_screening(ihe4, ibe8, yps, tk, rho, Fscreenaux)
            v(i) = v(i) * Fscreen * Fscreenaux

         case(109)

            ! n(2a,g)be9

            call calculate_screening(ihe4, ihe4, yps, tk, rho, Fscreen)
            v(i) = v(i) * Fscreen

         case default

            ! everything else

            if (k2(i) == 2) then

               ! catch reaction between idendical ions (e.g. 2(c12) -> something
               ! as opposed to c12 + c12 -> something)

               call calculate_screening(k1(i), k1(i), yps, tk, rho, Fscreen)

            else

               ! typical case

               call calculate_screening(k1(i), k3(i), yps, tk, rho, Fscreen)

            end if

            v(i) = v(i) * Fscreen

         end select
      end do

      tti2 = wallclocktime()

#ifdef PROFILE
      write (*,fmtprofile)'screen cpu_time/s = ' , tti2 - tti1
#endif


      ! turn off weak rates if requested
      if (.not. weak_rates) then
         do i = 1, nre
            if (.not. considerreaction(i)) cycle
            select case(ilabb(i))
            ! RJS 23/11/19 -- added i_bp
            case(i_ec,i_bm,i_bn,i_bp,i_ve,i_vp,i_vsp,i_vsn,i_vsa)
               v(i) = 1.e-99_r8
            end select
         end do
      end if

      ! turn off strong rates if requested
      if (.not. strong_rates) then
         do i = 1, nre
            if (.not. considerreaction(i)) cycle
            select case(ilabb(i))
            ! RJS 23/11/19 -- added i_bp
            case(i_ec,i_bm,i_bn,i_bp,i_ve,i_vp,i_vsp,i_vsn,i_vsa)
               continue
            case default
               v(i) = 1.e-99_r8
            end select
         end do
      end if

      ! use only weak rates, beta-delayed neutron/proton emission and alpha-decays if we're decaying
      ! TODO: should add p- and n- decays
      if (decay) then
         do i = 1, nre
            if (.not. considerreaction(i)) cycle
            select case(ilabb(i))
            ! RJS 23/11/19 -- added i_bp
            case(i_ec,i_bm,i_bn,i_bp,i_ba)
               cycle
            case default
               v(i) = 1.e-99_r8
            end select
         end do
      end if

      ! clip
      v = max(v,1.e-99_r8)

   end subroutine evaluate_all_rates

   subroutine apply_rate_factors()
         use physics_knobs, only: num_rate_factors, rate_index, rate_factor
         use reaction_info, only: rfac
         use rates, only: v
         implicit none
         integer(i4) :: i

         do i = 1, num_rate_factors
            if (rate_index(i) == -1) then
             cycle
            else
             v(rate_index(i)) = max(v(rate_index(i)) * rate_factor(i), 1.e-99_r8)
             if (abs(rfac(i)-1.) .gt. 1e-5) then
              stop 'you need serious help... you applied correction factor in both ppn_physics.input and in networksetup.txt??'
             end if
            end if 
         end do

         ! old rate factors from networksetup:
         v(:) = max(v(:) * rfac(:), 1.e-99_r8)
   end subroutine apply_rate_factors

end module evaluate_rates













