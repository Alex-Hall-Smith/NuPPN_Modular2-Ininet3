! ^_^
! This module provides routines for computing/interpolating rates that are not (yet) part of one of
! the standard reaction libraries that we use.
! 
! Current Limitations:
! -------------------
! 
! 1. Reactions can be inserted here in an overwrite mode, meaning that if you want to add a new
! reaction rate then we already have to have it in the network from some other source, and your rate
! will overwrite it, if activated.
!
! 2. Only reactions of the 15 types defined in reaction_info.F90 can be added. Things like 3a or
! c12+c12 or o16+o16 are currently handled by vital.F90. Sorry.

module other_nuc
   use utils
   use reaction_info
   use interp
   use akima
   use constants
   implicit none
   private
   
   ! ^_^ SWITCHES
   logical :: use_akima = .true. !< set true to use akima interpolation, otherwise linear

   integer(i4), public :: num_other
   integer(i4), allocatable, public :: &
      other_a(:), other_z(:), other_rtype(:)
   real(r8), public, allocatable :: other_rate(:)
   character(len=5), allocatable, public :: other_label(:)

   character(len=256), parameter :: &
      other_fn = "./other_reactions.txt", &
      other_fn_default = "../NPDATA/other_reactions.txt"

   ! ntrans for merging. hate this stuff
   integer(i4), public, allocatable :: other_nuc_ntrans(:,:,:,:)

   ! linear rate table interpolants
   type(interpolant) :: &
      w15a_ni63_intrp, &
      w17a_cu63_intrp, &
      p19a_cu65_intrp

   ! akima rate table interpolants
   type(akima_interpolant) :: &
      w15a_ni63_akima, &
      w17a_cu63_akima, &
      p19a_cu65_akima

   public :: other_nuc_init, calculate_other_nuc

contains

   subroutine other_nuc_init
      use communication
      implicit none
      integer(i4) :: i, other_fh
      character(len=256) :: fn
      logical :: user_file_exists

      ! read number of other reactions to include and their labels from disk
      if (master) then
         ! try to find user's version in run directory
         fn = other_fn_default
         inquire(file = trim(other_fn), exist = user_file_exists)
         if (user_file_exists) fn = other_fn
         open(file = fn, newunit = other_fh, status = "old", action = "read")
         read(other_fh,"(i4)") num_other
      end if

#ifndef PPN
      call broadcast(num_other)
#endif

      allocate( &
         other_a(num_other), &
         other_z(num_other), &
         other_rtype(num_other), &
         other_label(num_other), &
         other_rate(num_other) &
         )

      if (master) then
         do i = 1, num_other
            read(other_fh,"(A5)",end = 101) other_label(i)
         end do
         close(other_fh)
      end if

#ifndef PPN
      call broadcast_ch_arr(other_label)
#endif

      ! start adding reactions
      do i = 1, num_other
         select case(other_label(i))
         case("w15ng")
            call init_w15a_ni63(other_a(i), other_z(i), other_rtype(i), other_label(i))
         case("w17ng")
            call init_w17a_cu63(other_a(i), other_z(i), other_rtype(i), other_label(i))
         case("p19ng")
            call init_p19a_cu65(other_a(i), other_z(i), other_rtype(i), other_label(i))
         case default
            print *, "other_nuc: reaction '", other_label(i), "' in file 'other_reactions.txt'", &
              " does not have a corresponding init or eval routine"
            stop "check other_reactions.txt and other_nuc.F90" 
         end select
      end do

      ! allocate the stupid ntrans
      allocate( &
         other_nuc_ntrans( &
         minval(other_a):maxval(other_a), & ! dim 1
         minval(other_z):maxval(other_z), & ! dim 2
         1, &                               ! dim 3
         maxval(other_rtype) &              ! dim 4
         ) &
      )

      do i = 1, num_other
         other_nuc_ntrans(other_a(i),other_z(i),1,other_rtype(i)) = i
      end do

      return

101   stop "number of reactions in line 1 of other_reactions.txt is not correct"

   end subroutine other_nuc_init

   subroutine calculate_other_nuc(t9,rho,ye)
      implicit none
      integer(i4) :: i
      real(r8) :: t9, rho, ye
      
      do i = 1, num_other
         select case(other_label(i))
         case("w15ng")
            call eval_w15a_ni63(t9,rho,ye,other_rate(i))
         case("w17ng")
            call eval_w17a_cu63(t9,rho,ye,other_rate(i))
         case("p19ng")
            call eval_p19a_cu65(t9,rho,ye,other_rate(i))
         end select
      end do

   end subroutine calculate_other_nuc


   ! PHYSICAL REVIEW C 92, 045810
   ! 63Ni(n,γ) cross sections measured with DANCE
   ! Weigand et al. (2015)
   subroutine init_w15a_ni63(a,z,rtype,label)
      implicit none
      integer(i4) :: a, z, rtype
      character(len=5) :: label
      real(r8) :: &
         lrate(11), &
         EkeV(11) = [5,10,15,20,25,30,40,50,60,80,100], &
         aSEF(11) = [ 1.00_r8, 1.00_r8, 1.00_r8, 1.00_r8, 1.00_r8, 0.99_r8, 0.98_r8, 0.98_r8, &
         0.98_r8, 1.00_r8, 1.02_r8 ], &
         macs(11) = [ 220.62_r8, 121.87_r8, 92.58_r8, 78.33_r8, 69.33_r8, 62.68_r8, 52.71_r8, &
         45.21_r8, 39.29_r8, 30.59_r8, 24.56_r8 ]
      real(r8), parameter :: mu_g = amu*63._r8/64._r8

      ! convert MACS to a rate (cm^3 mol^-1 s^-1),
      ! e.g. Eq. 6 from Weigand et al. 2017
      lrate(:) = log10(macs*1e-27_r8*avogadro*aSEF*(2*EkeV(:)*kev2erg/mu_g)**HALF)

      a = 63; z = 28; rtype = i_ng; label = "w15ng"

      if (use_akima) then
         call create_akima_interpolant(EkeV, lrate, w15a_ni63_akima)
      else
         call create_1d_interpolant(11, EkeV, lrate, w15a_ni63_intrp)
      end if

   end subroutine init_w15a_ni63

   subroutine eval_w15a_ni63(t9,rho,ye,rate)
      implicit none
      real(r8), intent(in) :: t9, rho, ye
      real(r8) :: EkeV, lrate, rate

      EkeV = t9 * 1.e9_r8 * boltzmann / kev2erg

      if (use_akima) then
         call interp_akima(w15a_ni63_akima, EkeV, lrate)
      else
         call interpolate_1d(w15a_ni63_intrp, EkeV, lrate)
      end if

      rate = 10._r8**lrate * rho

   end subroutine eval_w15a_ni63


   ! NEXT RATE PLEASE.....

   ! PHYSICAL REVIEW C 95, 015808
   ! 63Cu(n,γ) cross section measured via 25 keV activation and time of flight
   ! Weigand et al. (2017)
   subroutine init_w17a_cu63(a,z,rtype,label)
      implicit none
      integer(i4) :: a, z, rtype
      character(len=5) :: label
      real(r8) :: &
         lrate(11), &
         EkeV(11) = [5,10,15,20,25,30,40,50,60,80,100], &
         macs(11) = [260.3_r8, 166.3_r8, 125.7_r8, 104.6_r8, 92.1_r8, 84.0_r8, &
            74.1_r8, 68.0_r8, 63.5_r8, 56.9_r8, 51.8_r8 ]
      real(r8), parameter :: mu_g = amu*63._r8/64._r8

      ! convert MACS to a rate (cm^3 mol^-1 s^-1),
      ! e.g. Eq. 6 from Weigand et al. 2017 (alpha_SEF = 1)
      lrate(:) = log10(macs*1e-27_r8*avogadro*(2*EkeV(:)*kev2erg/mu_g)**HALF)

      a = 63; z = 29; rtype = i_ng; label = "w17ng"

      if (use_akima) then
         call create_akima_interpolant(EkeV, lrate, w17a_cu63_akima)
      else
         call create_1d_interpolant(11, EkeV, lrate, w17a_cu63_intrp)
      end if

   end subroutine init_w17a_cu63

   subroutine eval_w17a_cu63(t9,rho,ye,rate)
      implicit none
      real(r8), intent(in) :: t9, rho, ye
      real(r8) :: EkeV, macs, lrate, rate

      EkeV = t9 * 1.e9_r8 * boltzmann / kev2erg

      if (use_akima) then
         call interp_akima(w17a_cu63_akima, EkeV, lrate)
      else
         call interpolate_1d(w17a_cu63_intrp, EkeV, lrate)
      end if

      rate = 10._r8**lrate * rho

   end subroutine eval_w17a_cu63


   ! In preparation - Prokop et al (2019) Cu65(n,g)
   subroutine init_p19a_cu65(a,z,rtype,label)
      implicit none
      integer(i4) :: a, z, rtype
      character(len=5) :: label
      real(r8) :: &
         lrate(11), &
         EkeV(11) = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100], &
         macs(11) = [ 1.32584e+02, 7.88296e+01, 5.81816e+01, 4.76608e+01, &
         4.13137e+01, 3.70363e+01, 3.15265e+01, 2.80125e+01, 2.54974e+01, &
         2.20207e+01, 1.96581e+01 ]
         ! older data again...
         !EkeV(11) = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100], &
         !   macs(11) = [ 1.34346e+02, 8.03893e+01, 5.97680e+01, 4.92270e+01, &
         !   4.28633e+01, 3.85857e+01, 3.31141e+01, 2.96650e+01, 2.72234e+01, &
         !   2.38929e+01, 2.16601e+01 ]
         ! data sent from Chris Prokop on 2018/11/08:
         !EkeV(12) = [5, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100], &
         !macs(12) = [ 1.66846e+02_r8, 9.54319e+01_r8, 5.28587e+01_r8, 4.46262e+01_r8, &
            !3.93776e+01_r8, 3.30554e+01_r8, 2.92195e+01_r8, 2.64850e+01_r8, &
            !2.43471e+01_r8, 2.25861e+01_r8, 2.10895e+01_r8, 1.97913e+01_r8 ]
         ! second (?) round of data from Chris Prokop:
         !EkeV(12) = [5, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100], &
         !macs(12) = [ 1.71649e+02_r8, 9.61047e+01_r8, 5.20037e+01_r8, 4.37388e+01_r8, &
            !3.85644e+01_r8, 3.25190e+01_r8, 2.90275e+01_r8, 2.66529e+01_r8, &
            !2.48668e+01_r8, 2.34381e+01_r8, 2.22497e+01_r8, 2.12339e+01_r8 ]
      real(r8), parameter :: mu_g = amu*65._r8/66._r8

       ! convert MACS to a rate (cm^3 mol^-1 s^-1),
       ! e.g. Eq. 6 from Weigand et al. 2017
       lrate(:) = log10(macs*1e-27_r8*avogadro*(2*EkeV(:)*kev2erg/mu_g)**HALF)
 
       a = 65; z = 29; rtype = i_ng; label = "p19ng"
 
       if (use_akima) then
          call create_akima_interpolant(EkeV, lrate, p19a_cu65_akima)
       else
          call create_1d_interpolant(11, EkeV, lrate, p19a_cu65_intrp)
      end if

   end subroutine init_p19a_cu65

   subroutine eval_p19a_cu65(t9,rho,ye,rate)
      implicit none
      real(r8), intent(in) :: t9, rho, ye
      real(r8) :: EkeV, lrate, rate
      real(r8), parameter :: mu_g = amu*65._r8/66._r8
      EkeV = t9 * 1.e9_r8 * boltzmann / kev2erg

      if (use_akima) then
         call interp_akima(p19a_cu65_akima, EkeV, lrate)
      else
         call interpolate_1d(p19a_cu65_intrp, EkeV, lrate)
      end if

      rate = 10._r8**lrate * rho

   end subroutine eval_p19a_cu65


end module other_nuc



