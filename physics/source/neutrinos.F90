! MODULE: neutrinos
!
!> @author Samuel Jones
!
!> Currently we only have the following capture and spallation reactions, although this could be
!> generalised to include spallation of heavier isotopes:
!> 
!> capture (charged current):
!> nu_e    + (a,z)   -> (a,z+1)   + e-
!> nubar_e + (a,z)   -> (a,z-1)   + e+ 
!>
!> spallation (neutral current):
!> nu_x + (a,z) -> (a,z-1)   + p  + nu_x'
!> nu_x + (a,z) -> (a-1,z)   + n  + nu_x'
!> nu_x + (a,z) -> (a-2,z-2) + a  + nu_x'

module neutrinos
   use utils
   use communication
   use physics_knobs
   use constants
   use nuc_data
   use akima
   implicit none

   integer(i4), parameter, public :: &
      i_nue  = 0, & !< electron neutrino capture
      i_nup  = 1, & !< antielectron neutrino capture
      i_nusp = 2, & !< spallation of proton
      i_nusn = 3, & !< spallation of neutron
      i_nusa = 4    !< spallation of alpha particle

   integer(i4), public, allocatable :: &
      nu_type(:), & !< type of neutrino reaction this is
      nu_a(:), &    !< target mass numbers
      nu_z(:), &    !< target proton numbers
      i_nu          !< location of start of neutrino reactions in network

   integer, public :: &
      num_packed_nu, &
      num_nu_reac !< number of neutrino reactions total

   character(len=5), public, allocatable :: &
      nu_par(:), &
      nu_dau(:)

   character(len=5), public :: &
      nu_name = "NTRNO"

   character(len=*), parameter :: &
      w90nc_fn = "../NPDATA/neutrinos/woosley1990_signeutral.txt", &
      w90nb_fn = "../NPDATA/neutrinos/woosley1990_br.txt", &
      w90cc_fn = "../NPDATA/neutrinos/woosley1990_sigcharged.txt"

   real(r8), allocatable :: &
      w90nc(:,:), & !< woosley et al 1990 neutral current cross sections vs temperature
      w90nb(:,:), & !< woosley et al 1990 neutral current branching ratios (spallation)
      w90cc(:,:), & !< woosley et al 1990 charged current cross sections vs temperature
      w90nT(:), &   !< woosley et al 1990 neutral current temperature coordinates (MeV)
      w90cT(:), &   !< woosley et al 1990 neutral current temperature coordinates (MeV)
      w90nv(:), &   !< woosley et al 1990 neutral current rates (once calculated)
      w90cv(:)      !< woosley et al 1990 charged current rates (once calculated)

   integer(i4), allocatable :: &
      w90na(:), & !< woosley et al 1990 target mass numbers (neutral)
      w90nz(:), & !< woosley et al 1990 target charge numbers (neutral)
      w90ca(:), & !< woosley et al 1990 target mass numbers (charged)
      w90cz(:)    !< woosley et al 1990 target proton numbers (charged)

   integer(i4), parameter :: &
      nTn = 5, & !< number of temperature points in woosley 1990 cross sections (neutral)
      nTc = 4    !< number of temperature points in woosley 1990 cross sections (charged)

   integer(i4) :: num_neutral, num_charged

   real(r8) :: &
      mf96v(2) !< neutrino capture rates on free nucleons from McLaughlin & Fuller 1996

   type(akima_interpolant), allocatable :: &
      w90nc_ai(:), & !< interpolant for neutral
      w90cc_ai(:)    !< interpolant for charged
      
   character(len=5), allocatable :: &
      w90npzis(:), & !< woosley et al 1990 neutral current target names
      w90cpzis(:)    !< woosley et al 1990 charged target names

   real(r8), public, allocatable :: &
      v_nu(:), &
      v_nu_packed(:)

   logical, public, allocatable :: &
      nu_mask(:)

   type, public :: neutrino
      real(r8) :: &
         lnue,  &   !< electron neutrino luminosity (from proto-NS)
         lnuebar, & !< anti-electron neutrino luminosity (from proto-NS)
         lnux,  &   !< "x"-neutrino (tau and mu) luminosity (from proto-NS)
         lnuxbar, & !< anti-"x" neutrino luminosity (from proto-NS)
         enue, &    !< average energy of the electron neutrino spectrum
         enuebar, & !< average energy of the antielectron neutrino spectrum
         enux, &    !< average energy of the "x" neutrino spectrum
         enuxbar, & !< average energy of the anti-"x" neutrino spectrum
         r, &       !< radius of material from proto-NS center
         tau_nu, &  !< e-folding time of neutrino luminosity
         t          !< time after shock breakout of the neutrinosphere
   end type neutrino

   logical, public :: nu_masks_exist = .false.

   public :: neutrinos_init, init_one_neutrino, calculate_neutrino_rates

contains

   subroutine init_one_neutrino(nu)
      implicit none
      type(neutrino) :: nu

      nu%lnue    = 1.e-99_r8
      nu%lnuebar = 1.e-99_r8
      nu%lnux    = 1.e-99_r8
      nu%lnuxbar = 1.e-99_r8
      nu%enue    = 1.e-99_r8
      nu%enuebar = 1.e-99_r8
      nu%enux    = 1.e-99_r8
      nu%enuxbar = 1.e-99_r8
      nu%r       = 1.e-99_r8
      nu%t       = 1.e-99_r8
      nu%tau_nu  = 1.e-99_r8
   end subroutine init_one_neutrino

   subroutine neutrinos_init(nu)
#ifdef PPN
      use frame_knobs, only: lnue, lnuebar, lnux, lnuxbar, &
         enue, enuebar, enux, enuxbar, nu_t, tau_nu, r
#endif
      use frame_knobs
      implicit none
      type(neutrino) :: nu
      integer(i4) :: i, cnt

      call init_one_neutrino(nu)

      if (.not. do_neutrinos) then
         return
      end if

      ! for the single-zone frame we can just take the neutrino properties as constant,
      ! until such time as we put in trajectories that have this information
#ifdef PPN
      nu%lnue    = lnue
      nu%lnuebar = lnuebar
      nu%lnux    = lnux
      nu%lnuxbar = lnuxbar
      nu%enue    = enue
      nu%enuebar = enuebar
      nu%enux    = enux
      nu%enuxbar = enuxbar
      nu%r       = r
      nu%t       = nu_t
      nu%tau_nu  = tau_nu
#endif

      call read_woosley_1990()

      ! account for the free nucleon reactions from McLaughlin & Fuller 1996
      ! in the total number of neutrino reactions (that's the +2)
      ! NB there are 2 charged reactions for each target (nu and nubar) and 3
      ! neutral current reactions (p, n and alpha exit channels)
      num_nu_reac = num_charged*2 + num_neutral*3 + 2

      allocate( &
         nu_par(num_nu_reac), &
         nu_dau(num_nu_reac), &
         nu_mask(num_nu_reac), &
         nu_type(num_nu_reac), &
         v_nu(num_nu_reac), &
         nu_a(num_nu_reac), &
         nu_z(num_nu_reac) &
         )

      ! set up capture and spallation reactions from woosley 1990
      do i = 1, num_charged
         nu_a(i) = w90ca(i)
         nu_z(i) = w90cz(i)
         nu_type(i) = i_nue
      end do
      do i = num_charged+1, num_charged*2
         nu_a(i) = w90ca(i)
         nu_z(i) = w90cz(i)
         nu_type(i) = i_nup
      end do
      cnt = num_charged*2
      do i = 1, num_neutral
         nu_a(cnt+1) = w90na(i); nu_z(cnt+1) = w90nz(i); nu_type(cnt+1) = i_nusn
         nu_a(cnt+2) = w90na(i); nu_z(cnt+2) = w90nz(i); nu_type(cnt+2) = i_nusp
         nu_a(cnt+3) = w90na(i); nu_z(cnt+3) = w90nz(i); nu_type(cnt+3) = i_nusa
         cnt = cnt + 3
      end do

      ! lastly, add the neutrino capture by free protons and neutrons
      nu_a(cnt+1) = 1; nu_z(cnt+1) = 1; nu_type(cnt+1) = i_nup
      nu_a(cnt+2) = 1; nu_z(cnt+2) = 0; nu_type(cnt+2) = i_nue

      v_nu(:) = 1.e-99_r8

   end subroutine neutrinos_init

   subroutine read_woosley_1990()
      implicit none
      integer(i4) :: i, w90nc_fh, w90nb_fh, w90cc_fh
      character(len=*), parameter :: &
         w90nc_fmt = "(2(i2,2x),a5,1x,5(es9.2))", &
         w90nb_fmt = "(2(i2,2x),a5,1x,4(es9.2))", &
         w90cc_fmt = "(2(i2,2x),a5,1x,4(es9.2))"

      if (master) then
         open(file=w90nc_fn, newunit=w90nc_fh, status="old", action="read")
         open(file=w90nb_fn, newunit=w90nb_fh, status="old", action="read")
         open(file=w90cc_fn, newunit=w90cc_fh, status="old", action="read")

         ! neutral current cross sections and branchings
         ! skip headers
         do i = 1, 8
            read(w90nc_fh,*)
            read(w90nb_fh,*)
         end do

         ! read number of targets
         read(w90nc_fh,*) num_neutral
         read(w90nb_fh,*)

         ! charged current cross sections
         ! skip headers
         do i = 1, 4
            read(w90cc_fh,*)
         end do

         ! read number of targets
         read(w90cc_fh,*) num_charged
      endif

      ! broadcast size and allocate data arrays
#ifndef PPN
      call broadcast(num_neutral); call broadcast(num_charged)
#endif

      allocate( &
         w90nc(num_neutral,nTn), &
         w90nb(num_neutral,4), &
         w90npzis(num_neutral), &
         w90na(num_neutral), &
         w90nz(num_neutral), &
         w90nT(nTn), &
         w90cc(num_charged*2,nTc), &
         w90ca(num_charged*2), &
         w90cz(num_charged*2), &
         w90cpzis(num_charged*2), &
         w90cT(nTc), &
         w90nv(num_neutral), &
         w90cv(num_charged*2), &  !< *2 for neutrino and anti-neutrino reactions on same target
         w90nc_ai(num_neutral), &
         w90cc_ai(num_charged*2) &
         )

      ! read data 
      if (master) then
         ! neutral current cross sections and branching ratios
         read(w90nc_fh,*) w90nT(:)
         do i = 1, num_neutral
            read(w90nc_fh,w90nc_fmt) w90nz(i), w90na(i), w90npzis(i), w90nc(i,:)
            read(w90nb_fh,w90nb_fmt) w90nz(i), w90na(i), w90npzis(i), w90nb(i,:)
         end do

         ! charged current 
         read(w90cc_fh,*) w90cT(:)
         read(w90cc_fh,*) 
         do i = 1, num_charged
            read(w90cc_fh,w90cc_fmt) w90cz(i), w90ca(i), w90cpzis(i), w90cc(i,:)
         end do
         read(w90cc_fh,*) 
         do i = num_charged+1, 2*num_charged
            read(w90cc_fh,w90cc_fmt) w90cz(i), w90ca(i), w90cpzis(i), w90cc(i,:)
         end do
      end if

#ifndef PPN
      call broadcast(w90nz); call broadcast(w90na); call broadcast(w90nc)
      call broadcast(w90nb); call broadcast(w90cz); call broadcast(w90ca)
      call broadcast(w90cc); call broadcast(w90nT); call broadcast(w90cT)
      call broadcast_ch_arr(w90npzis); call broadcast_ch_arr(w90cpzis)
#endif

      ! convert cross sections to cm^2
      w90nc = w90nc * 1e-42
      w90cc = w90cc * 1e-42

      ! set up neutral cross section interpolants
      do i = 1, num_neutral
         call create_akima_interpolant(w90nT,w90nc(i,:),w90nc_ai(i))
      end do

      ! set up charged cross sectoin interpolants
      do i = 1, num_charged*2
         call create_akima_interpolant(w90cT,w90cc(i,:),w90cc_ai(i))
      end do

      if (master) then
         close(w90nc_fh); close(w90nb_fh); close(w90cc_fh)
      endif

   end subroutine read_woosley_1990

   subroutine nu_create_masks()
      implicit none
      integer(i4) :: i

      nu_mask(:) = .false.
      do i = 1, num_nu_reac
         select case(nu_type(i))
         case(i_nue)
            if (ispe(epsi(nu_a(i),nu_z(i))) /= -1 .and. ispe(epsi(nu_a(i),nu_z(i)+1)) /= -1) then
               nu_mask(i) = .true.
            end if
         case(i_nup)
            if (ispe(epsi(nu_a(i),nu_z(i))) /= -1 .and. ispe(epsi(nu_a(i),nu_z(i)-1)) /= -1) then
               nu_mask(i) = .true.
            end if
         case(i_nusp)
            if (ispe(epsi(nu_a(i),nu_z(i))) /= -1 .and. ispe(epsi(nu_a(i)-1,nu_z(i)-1)) /= -1) then
               nu_mask(i) = .true.
            end if
         case(i_nusn)
            if (ispe(epsi(nu_a(i),nu_z(i))) /= -1 .and. ispe(epsi(nu_a(i)-1,nu_z(i))) /= -1) then
               nu_mask(i) = .true.
            end if
         case(i_nusa)
            if (ispe(epsi(nu_a(i),nu_z(i))) /= -1 .and. ispe(epsi(nu_a(i)-2,nu_z(i)-2)) /= -1) then
               nu_mask(i) = .true.
            end if
         end select
      end do

      num_packed_nu = count(nu_mask)

      allocate(v_nu_packed(num_packed_nu))
      nu_masks_exist = .true.
   end subroutine nu_create_masks

   ! evaluate neutrino reaction rates from cross sections and other information about the neutrinos
   subroutine w90_calculate_rates(nu)
      use utils, only: bsearch_r8
      implicit none
      type(neutrino) :: nu
      real(r8) :: &
         nc(num_neutral), cc(num_charged*2), & !< interpolated neutral and charged cross sections
         Tnue, Tnuebar, Tnux, Tnuxbar, & !< electron and "x" neutrino temperatures (MeV)
         fnue, fnuebar, fnux, fnuxbar, & !< neutrino fluxes
         prefac
      integer(i4) :: i

      Tnue = nu%enue / 3.15_r8; Tnuebar = nu%enuebar / 3.15_r8
      Tnux = nu%enux / 3.15_r8; Tnuxbar = nu%enuxbar / 3.15_r8

      ! neutrino fluxes (woosley et al. 1990 equation 5)
      prefac = 9.9e56_r8/(4*pi*nu%r**2) * 10._r8 * &
         (ONE/nu%tau_nu) * exp(-nu%t / nu%tau_nu)

      ! I think Woosley assumes that the energy is equally partitioned between all types of
      ! neutrinos, but we have luminosities as input so we can adjust for this here
      fnue    = nu%lnue    / (nu%lnue + nu%lnuebar + nu%lnux + nu%lnuxbar) * prefac / Tnue
      fnuebar = nu%lnuebar / (nu%lnue + nu%lnuebar + nu%lnux + nu%lnuxbar) * prefac / Tnuebar
      fnux    = nu%lnux    / (nu%lnue + nu%lnuebar + nu%lnux + nu%lnuxbar) * prefac / Tnue
      fnuxbar = nu%lnuxbar / (nu%lnue + nu%lnuebar + nu%lnux + nu%lnuxbar) * prefac / Tnuebar

      ! interpolate cross sections
      do i = 1, num_neutral
         call interp_akima(w90nc_ai(i), HALF*(Tnux+Tnuxbar), nc(i))
         w90nv(i) = nc(i) * (fnue + fnuebar + fnux + fnuxbar) ! rate (/s) not including branching
      end do
      do i = 1, num_charged
         call interp_akima(w90cc_ai(i), Tnue, cc(i))
         w90cv(i) = cc(i) * fnue
         call interp_akima(w90cc_ai(num_charged+i), Tnuebar, cc(num_charged+i))
         w90cv(num_charged+i) = cc(num_charged+i) * fnuebar ! rate (/s) not including branching
      end do

   end subroutine w90_calculate_rates

   subroutine mf96_calculate_rates(nu)
      implicit none
      type(neutrino) :: nu
      real(r8) :: &
         c1, c2, &      !< terms from eqs 14 and 15 in Mclaughlin & Fuller 1996
         Tnue, Tnuebar, & !< electron and antielectron neutrino temperatures (MeV)
         r72inv

      ! the radius appears as r_7^2, implying units for r_7 of 1e7 cm, however the
      ! units of the rates lambda_* are s^-1, implying that r_7 is unitless and is
      ! r_7 = r / (1e7 cm)
      r72inv = ONE / ((nu%r/1e7_r8)**2)

      ! estimate temperature of blackbody distribution to be (Burrows & Lattimer 1987 ApJL)
      Tnue    = nu%enue    / 3.15_r8
      Tnuebar = nu%enuebar / 3.15_r8

      ! constants from Mclaughlin & Fuller 1996 eqs 14 and 15
      c1 = ONE + (0.6283_r8/Tnue)   + (0.1292_r8/Tnue**2)
      c2 = ONE + (1.158_r8/Tnuebar) + (0.6_r8/Tnuebar**2) + (0.1778_r8/Tnuebar**3)

      ! nu_e    + n   -> p   + e- (Mclauglin & Fuller 1996 eq 12)
      mf96v(1) = 0.1945_r8 * (nu%lnue/1e51_r8) * Tnue * r72inv * c1

      ! nubar_e + p   -> n   + e+ (Mclauglin & Fuller 1996 eq 13)
      mf96v(2) = 0.2_r8 * exp(-1.804_r8/Tnue) * (nu%lnuebar/1e51_r8) * Tnuebar * r72inv * c2

   end subroutine mf96_calculate_rates

   !> Given a electron and antielectron neutrino luminosities lnu and lnubar from the proto-neutron
   !> star and the radius of the particle being exposed, we calculate the neutrino capture rates
   !> on free nucleons after McLaughlin & Fuller 1996 ApJ 472 440
   !> In PPN, the neutrino luminosities and the radius are given in the frame input deck
   subroutine calculate_neutrino_rates(t9, rho, ye, nu)
      implicit none
      integer :: i, k, cnt
      real(r8) :: t9, rho, ye
      type(neutrino) :: nu

      ! calculate reaction rates from Woosley et al. 1990
      call w90_calculate_rates(nu)

      ! construct neutrino reaction rate array from calculated rates (charged current)
      do i = 1, num_charged*2
         v_nu(i) = w90cv(i)
      end do
      cnt = num_charged*2
      ! and rates+branchings (for neutral current)
      do i = 1, num_neutral
         v_nu(cnt+1) = w90nv(i)*w90nb(i,1)
         v_nu(cnt+2) = w90nv(i)*w90nb(i,2)
         v_nu(cnt+3) = w90nv(i)*w90nb(i,3)
         cnt = cnt + 3
      end do

      ! now for neutrino capture by free nucleons from McLaughlin & Fuller 1996
      call mf96_calculate_rates(nu)

      v_nu(cnt+1) = mf96v(1)
      v_nu(cnt+2) = mf96v(2)

      v_nu = max(v_nu,1e-99_r8)

      k = 0
      do i = 1, num_nu_reac
         if (nu_mask(i)) then
            k = k + 1
            v_nu_packed(k) = v_nu(i)
         end if
      end do

   end subroutine calculate_neutrino_rates

end module neutrinos


