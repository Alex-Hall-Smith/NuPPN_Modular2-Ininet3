!> @author Samuel Jones
!>
!> KADoNiS
!> -------
!> This module provides routines to read and interpolate reaction rates from
!> Bao et al 2000, Atomic Data and Nuclear Data Tables v76, p70-154
!> and
!> KADoNiS (Kalsruhe Astrophysical Database of Nucleosynthesis in Stars).
!>
!> These rates are experimental (n,g) rates (or the theoretical rates, if there are no
!> info from experiments) according to
!>
!> * MACS = Maxwellian-Averaged Cross Sections
!> * SEFs = Stellar Enhancement Factors

module kadonis
   use utils, only: r8, bsearch_r8
   use constants
   use communication
   use akima
   use interp
   use physics_knobs, only: kadonis_interp
   private

   ! ^_^ switches
   ! Interpolate the rate (Na<ov>) insead of cross section (RECOMMENDED).
   ! Interpolating cross sections and then introducing energy dependence can introduce unphysical
   ! oscillations in the reaction rate.
   ! note that if choosing to interpolate MACS, currently only linear interpolation is supported
   logical, parameter :: interpolate_rate = .true.

   ! extrapolate cross section/rate to high temperatures using 1/v (1/sqrt(2*Ekev/reduced_mass))
   ! scaling
   logical, parameter :: extrapolate_rate = .true.

   ! ^_^ parameters
   integer, public :: &
      nkad            !< number of (n,g) captures in file which are not duplicates
   !< (could be automated or written into file)
   integer :: &
      nsef    !< number of SEF in file seftable.txt

   integer, parameter :: &
      ntem = 11  !< number of columns/energies where MACS and SEF are given, extended if extrapolating.

   ! ^_^ file handles, length of file content, file names
   integer ::  macs_fh, sef_fh
   integer, parameter :: lenmacs = 449, lensef = 356
   ! ^_^ amount to extend tables (==0 if not extrapolating)
   integer :: nextra

   character(len=256), parameter :: &
      macs_fn      = '../NPDATA/macs_tablerev5updateK03.txt', &
      sef_fn       = '../NPDATA/seftable_K03.txt'

   real(r8) :: &
      kadonis_time  !< for profiling

   ! ^_^ tabulated data
   real(r8), allocatable :: &
      ener(:), & !< energy (keV) coordinates of MACS and SEF tables
      macs(:,:), &   !< cross-sections (mb)
      sefs(:,:), &   !< SEFs
      rate(:,:), &   !< rates (cm^3 mol^-1 s^-1)
      mu_g(:), &     !< reduced masses
      sef(:)         !< interpolated stellar enhancement factors
   real(r8), allocatable, public :: &
      vkad(:)        !< rates (output: interpolated)

   integer, allocatable, public :: &
      zkad(:), &      !< proton number in macs file
      akad(:)         !< mass number in macs file

   !> for storing interpolant data so we can do linear/Akima interpolation in energy/temperature
   type(akima_interpolant), allocatable, dimension(:) :: &
      aic, & ! cross section
      ais    ! stellar enhancement factors
   type(interpolant), allocatable, dimension(:) :: &
      lic, & ! cross section
      lis    ! stellar enhancement factors

   public :: kadonis_time, kadonis_init, kadonis_interpolate_rates

contains


   !---
   !> Initialization of kadonis module
   !---
   subroutine kadonis_init()
      implicit none
      integer(i4) :: i

      if (extrapolate_rate) then
         nextra = 13 ! 9 above up to 1 MeV, 4 below down to 1 keV
      else
         nextra = 0
      end if

      call kadonis_read_data()

      ! calculate reduced masses
      allocate(mu_g(nkad))
      mu_g(:) = amu * real(akad,r8) / real((akad+1),r8)

      ! extrapolate the cross sections, if requested
      call kadonis_extrapolate_macs()

      ! We now have MACS and SEF data for all of the reactions frrom kadonis at 11 energy points.
      ! Interpolating the rate is more physical than interpolating the cross section and adding
      ! an additional energy dependence when it is converted to a rate, so here we convert all the
      ! cross section data points to rates

      select case(kadonis_interp)
      case(1)
         allocate(lic(nkad))
         if (.not. interpolate_rate) allocate(lis(nkad))
      case(3)
         allocate(aic(nkad))
         if (.not. interpolate_rate) allocate(ais(nkad))
      end select

      do i = 1, nkad
         if (interpolate_rate) then
            ! rate = <sigma> * N_A * alpha_SEF * sqrt(2E/m)
            ! see, e.g. Bao 2000 or Weigand et al. 2017
            rate(i,:) = macs(i,:)*avogadro*sefs(i,:)*(2*ener(:)*kev2erg/mu_g(i))**HALF
         else
            rate(i,:) = log10(macs(i,:))
         end if

         select case(kadonis_interp)
         case(1) ! linear
            call create_1d_interpolant(ntem+nextra, ener, rate(i,:), lic(i))
            if (.not. interpolate_rate) call create_1d_interpolant(ntem+nextra, ener, sefs(i,:), lis(i))
         case(3) ! akima
            call create_akima_interpolant(ener, rate(i,:), aic(i))
            if (.not. interpolate_rate) call create_akima_interpolant(ener, sefs(i,:), ais(i))
         case default
            stop "kadonis: invalid kadonis_interp option in ppn_physics.input"
         end select
      end do

   end subroutine kadonis_init


   !---
   !> Interpolate the KADoNiS rates at a given energy
   !---
   subroutine kadonis_interpolate_rates(t9,rho)
         implicit none
         integer  :: i, j, k, khi, klo
         real(r8) :: Ekev_in, Ekev, Ehi, Elo, DE, &
            baosef(nkad)                  !< SEFs
         real(r8), intent(in) :: t9, rho

         Ekev_in = t9 * 1e9_r8 * (boltzmann/kev2erg)

         ! ^_^ clip energy to table
         Ekev = min(max(Ekev_in, ener(1)), ener(ntem+nextra))

         ! interpolate rates and multiply with density
         select case(kadonis_interp)
         case(1) ! linear
            do i = 1, nkad
               call interpolate_1d(lic(i), Ekev, vkad(i))
               if (.not. interpolate_rate) call interpolate_1d(lis(i), Ekev, sef(i))
            end do
         case(3) ! akima
            do i = 1, nkad
               call interp_akima(aic(i), Ekev, vkad(i))
               if (.not. interpolate_rate) call interp_akima(ais(i), Ekev, sef(i))
            end do
         case default
            stop "kadonis: invalid kadonis_interp option in ppn_physics.input"
         end select

         if (.not. interpolate_rate) then
            ! sam's calculation:
            vkad(:) = 10._r8**vkad(:) * avogadro*sef(:)*(2*Ekev*kev2erg/mu_g)**HALF
            ! old ppn version of this calculation:
            !vkad(:) = vkad(:)/1.e-27 * baosef(:) * sqrt(Ekev_in/30._r8) * (1.442e5_r8* ((dble(akad)+1)/dble(akad)) ** HALF )
         end if

         vkad(:) = vkad(:) * rho

   end subroutine kadonis_interpolate_rates



   !---
   !> Reads the KADoNiS data
   !---
   subroutine kadonis_read_data()
         implicit none
         character(len=2) :: el
         integer(i4) :: i, j, cnt, theo, &
            amacs(lenmacs), zmacs(lenmacs), asef(lensef), zsef(lensef)
         real(r8) :: sef_from_file(lensef,ntem)
         logical :: mask_macs(lenmacs)

         if (master) then
            open(newunit = macs_fh, file = macs_fn, status = 'old')

            ! first just read the (Z,A) for all isotopes in each file and make a mask that is .true.
            ! only for the rates that we want to use, which are the first-appearing in case of
            ! duplicates
            mask_macs(:) = .true.
            read(macs_fh,*) !< skip coordinates
            do i = 1, lenmacs
               read(macs_fh,*) zmacs(i), amacs(i)

               if (i == 1) cycle

               if (zmacs(i) == zmacs(i-1) .and. amacs(i) == amacs(i-1)) then
                  mask_macs(i) = .false.
               end if
            end do

            nkad = count(mask_macs)
         end if

         ! now we have the real dimension of the kadonis rates we will be using, so
         ! broadcast and allocate
#ifndef PPN
         call broadcast(nkad)
#endif

         allocate( &
            macs(nkad,ntem), &  !< cross-sections
            sefs(nkad,ntem), &  !< SEFs
            rate(nkad,ntem), &  !< SEFs
            ener(ntem),      &  !< energies
            vkad(nkad),      &  !< rates for output
            sef(nkad),       &  !< interpolated SEFs
            zkad(nkad),      &  !< charge numbers
            akad(nkad)       &  !< mass numbers
            )

         ! now let's have the master rewind the MACS file read only the relevant data
         if (master) then
            rewind(macs_fh)

            ! energy coordinates
            read (macs_fh,*) (ener(i),i=1,ntem)

            ! MACS data
            cnt = 0
            do i = 1, lenmacs
               if (.not. mask_macs(i)) then
                  read(macs_fh,*)
               else
                  cnt = cnt + 1
                  read(macs_fh,*) zkad(cnt), akad(cnt), el, (macs(cnt,j),j=1,ntem), theo
               end if
            end do

            close(macs_fh)
         end if

#ifndef PPN
         call broadcast(ener); call broadcast(zkad); call broadcast(akad); call broadcast(macs)
#endif

         ! now read the SEF data, which does not apparently line up with the isotopes in the MACS file, so
         ! we need to take only the SEFs that we have MACS for, and assume all SEFs we don't have
         ! (if any) are unity across the energy range.

         sef_from_file(:,:) = ONE

         if (master) then
            open(newunit = sef_fh,  file = sef_fn,  status = 'old')
            read (sef_fh,*) ! skip coordinates
            do i = 1, lensef
               read (sef_fh,*) zsef(i), asef(i), el, (sef_from_file(i,j),j=1,ntem)
            end do
            close(sef_fh)

            ! search for SEFs that match MACS data and save them in the SEF array
            do i = 1, nkad
               do j = 1, lensef
                  if (zkad(i) == zsef(j) .and. akad(i) == asef(j)) then
                     ! this the right SEF
                     sefs(i,1:ntem) = sef_from_file(j,1:ntem)
                     exit
                  end if
               end do
            end do
         end if

#ifndef PPN
         call broadcast(sefs)
#endif

         ! finally, convert cross sections to cm^2 from mb
         macs(:,:) = macs(:,:) * 1.e-27_r8

   end subroutine kadonis_read_data


   subroutine kadonis_extrapolate_macs
      implicit none
      integer :: i, j
      real(r8), allocatable :: macs_old(:,:), sefs_old(:,:), ener_old(:)
      real(r8) :: gradl, gradh

      if (.not. extrapolate_rate) return

      ! temporary storage for macs, sefs and energies from file
      allocate(macs_old(nkad,ntem), sefs_old(nkad,ntem), ener_old(ntem))
      macs_old = macs; sefs_old = sefs; ener_old = ener

      ! re-allocate macs sefs and ener for new size (including extrapolated
      ! data)
      deallocate(macs,sefs,ener,rate)
      allocate(macs(nkad,ntem+nextra),sefs(nkad,ntem+nextra),ener(ntem+nextra),rate(nkad,ntem+nextra))

      ! set up new energy array
      ener(5:ntem+4) = ener_old(:)
      ener(1:4) = [1, 2, 3, 4]
      ener(ntem+5:ntem+nextra) = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

      ! extrapolate up to 1 MeV (~10 GK) if requested
      do i = 1, nkad ! loop over reactions

         ! preserve macs and sefs from file for intermediate energies
         macs(i,5:ntem+4) = macs_old(i,:)
         sefs(i,5:ntem+4) = sefs_old(i,:)

         ! low side - gradients for log-linear extrapolation
         if (macs(i,5) == ZERO) then
            macs(i,1:4) = ZERO
         else
            gradl = ( log10(macs(i,6)) - log10(macs(i,5)) ) / ( ener(6) - ener(5) )
            do j = 1, 4
               macs(i,j) = 10._r8 ** (log10(macs(i,5)) - (ener(5) - ener(j))*gradl)
            end do
         end if
         sefs(i,1:4) = sefs(i,5)

         ! high side - gradients for log-linear extrapolation
         if (macs(i,ntem+4) == ZERO) then
            macs(i,ntem+5:) = ZERO
         else
            gradh = ( log10(macs(i,ntem+4)) - log10(macs(i,ntem+3)) ) / ( ener(ntem+4) - ener(ntem+3) )
            do j = ntem+5, ntem+nextra
               macs(i,j) = 10._r8 ** (log10(macs(i,ntem+4)) + (ener(j) - ener(ntem+4))*gradh)
            end do
         end if
         sefs(i,ntem+5:) = sefs(i,5)

      end do

      deallocate(macs_old,sefs_old,ener_old)

   end subroutine kadonis_extrapolate_macs



end module kadonis


