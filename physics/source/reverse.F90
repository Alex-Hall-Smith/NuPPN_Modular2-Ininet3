!> @Author Samuel Jones, Heiko Möller
!>
!> This module calculates the rate of reverse reactions using the principle of detailed
!> balance
module reverse
   use nse_swj, only: ppn_g, ppn_q, ppn_m, ppn2reaclib
   use reaclib, only: reaclib_interpolate_partition_functions, partinter
   use utils, only: r8, i4
   use array_sizes, only: nsp, nre
   use communication
   implicit none
   private

   integer(i4) :: irev(nre)
   logical :: this_is_a_reverse_reaction(nre)

   !> coulomb corrections need not be included, because we screen
   !> the reverse rate after the fact
   logical :: use_coulomb_corrections = .false.

   public calculate_reverse_rates, reverse_init

contains

   !> create a map of reaction to reverse reaction in the ppn network
   subroutine reverse_init()
         use rates, only: k1, k2, k3, k4, k5, k6, k7, k8
     ! RJS 23/11/19 -- added i_bp
         use reaction_info, only: ilabb, i_ec, i_bm, i_bn, i_bp, lab, i_ve, i_vp, i_vsp, i_vsn, i_vsa
         use nuc_data, only: considerreaction, isomeric_state, ispe
         use reaclib, only: is_reverse, reaclib_ntrans, aan, zzn
         implicit none
         integer(i4) :: i, j, ri(8), rj(8), isomer(8), iamgamma, ridx
         real(r8) :: anetw(nsp), znetw(nsp)
         common/cnetw/anetw,znetw

         if (master) print *, "init reverse rates ... "

         irev = -1
         this_is_a_reverse_reaction = .false.
         iamgamma = ispe("OOOOO")

         ! identify likely reverse reactions from reaclib
         do i = 1, nre
         !do i = 1, 2000
            if (.not. considerreaction(i)) cycle
! RJS 23/11/19 -- changed to 16 as new reaction type added
            if (ilabb(i) > 16) cycle !< non-standard reaction type, see reaction_info

            ! is this reaction a reverse reaction? - query reaclib
            ridx = reaclib_ntrans(int(anetw(k1(i))),int(znetw(k1(i))),1,ilabb(i))
            if (ridx == 0) cycle !< this reaction isn't in reaclib
            this_is_a_reverse_reaction(i) = is_reverse(ridx)
         end do

         do i = 1, nre
            if (.not. considerreaction(i)) cycle
            ! ignore weak rates and beta-delayed neutron/proton emission
            if (any(ilabb(i) == [i_ec,i_bm,i_bn,i_bp,i_ve,i_vp,i_vsp,i_vsn,i_vsa])) cycle

            ! ignore reactions with identical products: these are done in VITAL
            ! in an inconsistent way that breaks detailed balance
            if (k8(i) > 1) cycle

            ! ignore reaction if it is a reverse reaction - we don't need to know where
            ! *its* reverse is - that's doubling up the work un-necessarily
            if (this_is_a_reverse_reaction(i)) cycle
            
            do j = 1, nre
               if (.not. considerreaction(j)) cycle

               ri(:) = [k1(i),k2(i),k3(i),k4(i),k5(i),k6(i),k7(i),k8(i)]
               rj(:) = [k7(j),k8(j),k5(j),k6(j),k3(j),k4(j),k1(j),k2(j)]
               if ( all(ri == rj) ) then
                  ! the reverse reaction of i is j
                  ! make sure none of the reactants or products are isomers
                  isomer(:) = 1 ! init to ground state
                  if (k1(i) /= iamgamma) isomer(1) = isomeric_state(k1(i))
                  if (k3(i) /= iamgamma) isomer(2) = isomeric_state(k3(i))
                  if (k5(i) /= iamgamma) isomer(3) = isomeric_state(k5(i))
                  if (k7(i) /= iamgamma) isomer(4) = isomeric_state(k7(i))
                  if (k1(j) /= iamgamma) isomer(5) = isomeric_state(k1(j))
                  if (k3(j) /= iamgamma) isomer(6) = isomeric_state(k3(j))
                  if (k5(j) /= iamgamma) isomer(7) = isomeric_state(k5(j))
                  if (k7(j) /= iamgamma) isomer(8) = isomeric_state(k7(j))
                  if (all(isomer == 1)) then
                     ! the reaction doesn't involve isomers in any way

                     ! now make sure that the reverse itself is flagged as a reverse by
                     ! reaclib
                     if (this_is_a_reverse_reaction(j)) then
                        ! this is indeed considered a reverse reaction.
                        ! save this reverse so we can replace it with detailed balance later
                        irev(i) = j
                        ! set label for the same reference as fwd rate
                        lab(j) = "RVRSE"
                     end if
                  end if
                  ! exit nested loop since we found the reverse reaction
                  exit
               end if
            end do
         end do
   end subroutine reverse_init

   !> calculate reverse rates using detailed balance
   subroutine calculate_reverse_rates(rho, t9, ye)
      use nse_swj, only: ppn2nse, g_tot_i, mu_i_coul, coulomb_corrections, &
         partition_functions
      use rates, only: v, v_rev
      use reaction_info
      use constants, only: ZERO, ONE
      use nuc_data, only: zis, considerreaction
      use reaclib, only: map_ppn_to_reaclib, locate_vbranch_in_reaclib_rv, rflag2
      implicit none
      real(r8), intent(in) :: rho, t9, ye
      real(r8) :: G(nsp), muC(nsp), fac
      integer(i4) :: i, idx, ierr

      ! get partition functions from nse module (ultimately from reaclib for now)
      call partition_functions(t9)

      ! get coulomb corrections to the ion chemical potentials from the nse_swj module
      if (use_coulomb_corrections) call coulomb_corrections(t9, rho, ye)

      ! save them
      G = ONE ! initialize to one
      muC = ZERO ! initialize to zero
      do i = 1, nsp
         idx = ppn2nse(i)
         if (idx /= 0) then
            G(i) = g_tot_i(idx) ! these include the spin factors
            if (use_coulomb_corrections) muC(i) = mu_i_coul(idx)
         end if
      end do

      ! calculate reverse rates for all reverse reactions that have been identified as
      ! such from reaclib during reverse_init. do this using detailed balance
      do i = 1, nre

         if (.not. considerreaction(i)) cycle

         if (irev(i) /= -1) then

            ! this reaction's reverse is in the network, and reaclib has decreed that is
            ! is a true reverse rate (i.e. it should be computed from its corresponding
            ! forward rate using detailed balance).

            ! for very low forward rates (i.e. forward rate was clipped to ~1e-99
            ! to prevent numerical underflow or zeros), just use the regular
            ! reverse rate from whatever compilation it usually comes from, so
            ! that we don't actually use 1e-99 as the true forward rate. This
            ! has been shown to produce spurious abundances. One example was
            ! production of ho165 during the main sequence of a massive star.
            
            ! hence, only replace reverse rate with our version if fwd rate is > some small value, say 1e-97
            if (v(i) > 1.e-97_r8) then

               ! calculate the reverse reaction rate using detailed balance
               call reverse_rate(i, rho, t9, fac, G, muC, ierr)

               if (ierr /= 0) cycle

               ! replace its reverse rate with the one computed using detailed balance

               v_rev(i) = max(v(i) * fac, 1e-99_r8) ! reverse rate
               v(irev(i)) = v_rev(i)
            end if

         end if

      end do

   end subroutine calculate_reverse_rates

   !> calculate the rate of the reverse reaction to the reaction at index idx
   !> in the ppn network using detailed balance
   subroutine reverse_rate(idx, rho, t9, fac, G, muC, ierr)
         use constants, only: boltzmann, ONE, erg2mev, ZERO, TWO, pi, hbar_mev, &
               avogadro, planck_mev, planck_erg, clight2
         use rates, only: v, k1, k2, k3, k4, k5, k6, k7, k8
         use nuc_data, only: zis
         use reaction_info
         implicit none

         integer(i4) :: &
               i, &   !< counter
               a1, &  !< index of reactant 1 in ppn network
               a2, &  !< index of reactant 2 in ppn network
               b1, &  !< index of product 1 in ppn network
               b2, &  !< index of product 2 in ppn network
               na1, & !< number of a1 reactants
               na2, & !< number of a2 reactants
               nb1, & !< number of b1 reactants
               nb2, & !< number of b2 reactants
               ierr, &
               idx    !< index of reaction in ppn network

         real(r8), intent(in) :: &
               muC(nsp), & !< coulomb correction to ion chemical potential at (t9,rho,ye)
               G(nsp), &   !< partition functions at this temperature t9
               t9,   &     !< temperature in GK
               rho         !< density (cgs)

         real(r8) :: &
               kTMeV,&       ! < boltzmann * temperature
               kTerg,&       ! < boltzmann * temperature
               G_a1, &       ! < partition function of reactant 1
               G_a2, &       ! < partition function of reactant 2
               G_b1, &       ! < partition function of product 1
               G_b2, &       ! < partition function of product 2
               m_a1, &       ! < nuclear mass of reactant 1
               m_a2, &       ! < nuclear mass of reactant 2
               m_b1, &       ! < nuclear mass of product 1
               m_b2, &       ! < nuclear mass of product 2
               mu_a1, &      ! < nuclear mass of reactant 1
               mu_a2, &      ! < nuclear mass of reactant 2
               mu_b1, &      ! < nuclear mass of product 1
               mu_b2, &      ! < nuclear mass of product 2
               Q,    &       ! < reaction Q value
               rhofac, &     ! < for dividing rho out to give bare macs from rate
               fac           ! < reverse rate factor

         logical :: &
               use_partition_functions = .true.

         real(r8) :: &
               x, pfac, efac, mfac, mufac

         ierr = 0

         kTerg = boltzmann * t9 * 1.e9_r8
         kTMeV = kTerg * erg2mev

         ! indices of reactants and products in ppn network
         a1  = k1(idx) ; a2  = k3(idx) ; b1  = k5(idx) ; b2  = k7(idx)
         ! number of each reactant and product
         na1 = k2(idx) ; na2 = k4(idx) ; nb1 = k6(idx) ; nb2 = k8(idx)

         ! reactants
         
         G_a1 = G(a1)     ! Partition function (including spin factor)
         m_a1 = ppn_m(a1) ! masses
         mu_a1 = muC(a1)  ! coulomb corrections to ion chemical potentials (MeV)

         ! account for reactions where reactants are identical (k2 = 2) or not (k2 = 1)

         select case(na1)
         case(1)
            G_a2 = G(a2)
            m_a2 = ppn_m(a2)
            mu_a2 = muC(a2)
         case(2)
            G_a2 = G_a1
            m_a2 = m_a1
            mu_a2 = mu_a1
         case default
            ierr = 1
            return
         end select

         ! products

         G_b2  = G(b2)
         mu_b2 = muC(b2)
         m_b2 = ppn_m(b2)

         ! account for case where products are identical (k8 = 2) or not (k8 = 1)

         select case(nb2)
         case(1)
            ! account for either no further products (nb1 = 0) or a 2nd product (nb1 = 1)
            select case(nb1)
            case(1)
               G_b1  = G(b1)
               mu_b1 = muC(b1)
               m_b1 = ppn_m(b1)
            case(0)
               G_b1  = ZERO
               m_b1  = ZERO
               mu_b1 = ZERO
            end select
         case(2)
            G_b1  = G_b2
            m_b1  = m_b2
            mu_b1 = mu_b2
         case default
            ierr = 1
            return
         end select

         ! Calder et al. 2007 eq. A6 and A7
         select case(ilabb(idx))
         case(i_ng,i_pg,i_ag)
            ! partition function factor (including spin factors)
            pfac = G_a1 * G_a2 / G_b2

            ! mass/energy factor
            mfac = (m_a1 * m_a2 / m_b2 * kTerg * TWO * pi / planck_erg / planck_erg) ** 1.5_r8

            ! Q-value (of reverse rate)
            Q = -ppn_q(a1) - ppn_q(a2) + ppn_q(b2)

            ! density needs to be divided out for photodisintegrations
            rhofac = ONE / rho / avogadro

            ! coulomb factor
            mufac = exp((-mu_a1 - mu_a2 + mu_b2) / kTMeV)
         case default
            ! partition function factor (including spin factors)
            pfac = G_a1 * G_a2 / G_b1 / G_b2

            ! reduced mass factor
            mfac = (m_a1 * m_a2 / m_b1 / m_b2) ** 1.5_r8

            ! Q-value (of reverse rate)
            Q = -ppn_q(a1) - ppn_q(a2) + ppn_q(b1) + ppn_q(b2)

            rhofac = ONE

            ! coulomb factor
            mufac = exp((-mu_a1 - mu_a2 + mu_b2 + mu_b1) / kTMeV)
         end select

         ! Q-value and boltzmann factor, making sure we don't get infinities
         x = min(-Q / kTMeV, 200._r8)
         efac = exp(x)

         fac = mfac * efac * rhofac
         if (use_partition_functions) fac = fac * pfac
         if (use_coulomb_corrections) fac = fac * mufac

         fac = min(max(fac,1e-99_r8), 1e99_r8)

   end subroutine reverse_rate

end module reverse
