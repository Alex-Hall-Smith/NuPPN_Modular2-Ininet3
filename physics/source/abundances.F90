! MODULE: abundances
!
!> @author S. W. Jones and company
!
! DESCRIPTION:
!> @brief Interface to sets of abundnaces.
!> @details
!> Provides routines to load isotopic abundace data in various formats that have been used in the NuGrid collaboration and in the
!> broader community over the years. It also provides routines to initialise simulations with these initial abundances. The choice
!> of abundances to load is controlled by the parameter \p iabuini

!     iabuini = 1: read solar distir from gnetw.dat  
!     iabuini = 2: read free-formatted isotope name followed by abundance
!     iabuini = 3: hardwire some given abundances  
!     iabuini = 4: ye is hardwired  
!     iabuini = 5: read iso_massf* stype PPN file
!     iabuini = 11: initialise from iniab* style files (e.g. mppnp/USEEPP)  
!     iabuini = 12: read abundances from file massf_*.dat, produced by abund_at_masscoorinate  
!     iabuini = 100: read abundance from filename with the iniabu method in the utils python package  
!     iabuini = 101: for tppnp
!     iabuini = 102: plain isotope name and mass fraction

!> @param iabuini choice of how to initialize isotopic abundances. details in source

module abundances
   use utils
   use constants
   use frame_knobs, only: iabuini, ini_filename, ye_initial
   use nse_torch, only: testnse
   use nse_swj, only: compute_nse
#ifndef PPN
   use communication
#endif
   implicit none
   private

   ! ^_^ loaded abundance data:
   integer, public :: num_iso_in             !< number of isotopes for which abundances are loaded
   character(len=5), public :: iso_in(10000) !< names of isotopes we are loading
   real(r8), public :: x_in(10000)           !< mass fractions of isotopes loaded

   ! ^_^ file names and formats
   character(len=*), parameter :: &
         torino_fn  = '../NPDATA/gnetw.tab', &
         torino_fmt = '(a1,1x,a5,2(1x,e13.3))', &
         frisch_fmt = '(4x,a5,9x,d16.10)', &
         nugrid_fmt = '(5x,a5,4x,d17.10)', &
         ppn_fmt    = '(24x,ES11.5,2x,a5)', &
         free_fmt   = '(a5,es23.15)'

   integer, public :: external_file_handle !< handle of open file from external source (e.g. tppnp)

   public set_initial_abundances

contains

   !---
   !  DESCRIPTION:
   !> @brief   loads a set of abundances and sets them into the state
   !
   !> @param   y[out]      abundance array (mass fractions)
   !#> @var     id_network  position of isotopes in network
   !#> @var     id_loaded   porision in network of isotopes loaded from file
   !---
   subroutine set_initial_abundances( y )
         use array_sizes
         use nuc_data
         integer :: i, j, k, id_network(nsp), id_loaded(10000)
         real(r8) :: y(nsp)

         y = 1.e-99_r8
         call load_abundances

         ! hash the integers corresponding to species names for faster linear searching
         do i = 1, nsp
            id_network(i) = ispe(zis(i))
         end do
         do i = 1, num_iso_in
            ! id will be -1 if not in our network
            id_loaded(i)  = ispe(iso_in(i))
         end do

         j = num_iso_in
         do k = 1, nsp
            if (.not. considerisotope(k)) cycle
            do i = 1, j
               if (id_loaded(i) == id_network(k)) then
                  ! set abundance and pop out the one we have used
                  y         ( k ) = x_in      ( i )
                  iso_in    ( i ) = iso_in    ( j )
                  id_loaded ( i ) = id_loaded ( j )
                  x_in      ( i ) = x_in      ( j )
                  j               = j - 1
                  exit
               end if
            end do
         end do

         y = max(y, 1.e-99_r8)
   end subroutine set_initial_abundances


   !---
   !  DESCRIPTION:
   !> @brief reads abundances in some format or another into memory
   !---
   subroutine load_abundances
#if defined(PPN) || defined(MPPNP)
         print *, ' initialising abundances:'
         print *, '     iabuini = ', iabuini
         print *, '     ini_filename = ', ini_filename
#endif
         select case(iabuini)
         case(1)
            ! this was iabuini = 2 in mppnp ... ack
            call load_torino_trans_iron
         case(2)
            call load_iso_abu_free
         case(3)
            call load_by_hand
         case(4)
            call load_fix_ye_nse
         case(5)
            call load_ppn
         case(11)
            call load_urs_frischknecht_xin
         case(12)
            call load_abu_from_masscoordinate
         case(100)
            call load_nugridpy_utils_xin
         case(101)
            call load_tppnp
         case(102)
            call load_plain_ascii
         case default
            stop 'iabuini choice not recognised'
         end select
   end subroutine load_abundances


   subroutine load_iso_abu_free
      integer :: fh, cnt

      fh = get_unused_file_handle()
      open ( fh, file = ini_filename )

      cnt = 0
      do
         cnt = cnt + 1
         read ( fh, free_fmt, end = 1945 ) iso_in(cnt), x_in(cnt)
      end do
      1945 num_iso_in = cnt - 1
      close(fh)
   end subroutine load_iso_abu_free

   !---
   !  DESCRIPTION:
   !> @brief get abundances from a tppnp input file
   !---
   subroutine load_tppnp
         ! these are actually read in in the read_trajectory subroutine of the traj_data module
         continue
   end subroutine load_tppnp


   !---
   !  DESCRIPTION:
   !> @brief read abundances from a torino-format ascii file
   !---
   subroutine load_torino_trans_iron
         integer :: fh, cnt
         character(len=5) :: cdummy
         real(r8) :: fdummy

         fh = get_unused_file_handle()
         open ( fh, file = torino_fn )
         cnt = 1
         do while (.true.)
            read(fh,'(a1)', end=9343) cdummy
            if (cdummy == 'D') then
               backspace(fh)
               read( fh, torino_fmt) cdummy, iso_in(cnt), fdummy, x_in(cnt)
               cnt = cnt + 1
            end if
         end do
         9343 num_iso_in = cnt - 1
         close(fh)
   end subroutine load_torino_trans_iron


   !---
   !  DESCRIPTION:
   !> @brief abundances of species are hard-wired in this subroutine
   !---
   subroutine load_by_hand
         num_iso_in = 7

         iso_in ( 1 )  = 'HE  4'
         iso_in ( 2 )  = 'C  12'
         iso_in ( 2 )  = 'C  13'
         iso_in ( 4 )  = 'O  16'
         iso_in ( 5 )  = 'MG 25'
         iso_in ( 6 )  = 'NE 22'
         iso_in ( 7 )  = 'FE 56'

         x_in   ( 1 )  = 4.550e-1_r8
         x_in   ( 2 )  = 3.366e-1_r8
         x_in   ( 3 )  = 5.734e-2_r8
         x_in   ( 4 )  = 1.449e-1_r8
         x_in   ( 5 )  = 2.955e-4_r8
         x_in   ( 6 )  = 1.358e-2_r8
         x_in   ( 7 )  = 1.000e-3_r8
   end subroutine load_by_hand


   !---
   !  DESCRIPTION:
   !> @brief electron fraction \f$ Y_\mathrm{e} \f$ is hardwired into the \f$ ^{56}\mathrm{Ni} \f$ and neutron abundances
   !---
   subroutine load_fix_ye_nse
         real(r8) :: ye

         ye = ye_initial

         num_iso_in = 2
         iso_in ( 1 )  = 'NI 56'
         iso_in ( 2 )  = 'NEUT '
         x_in   ( 1 )  = ye / ( 28._r8 / 56._r8 ) 
         x_in   ( 2 )  = ONE - x_in(1) 

         !call compute_nse( t9, rho, ye, yps, mu_p, mu_n, iterations, ierr )
         
         write(*,*) 'IABUINI = 4: composition set to NSE for ye = ', ye_initial

   end subroutine load_fix_ye_nse

   !---
   !  DESCRIPTION:
   !> @brief ppn-style iso_massf* file is read
   !---
   subroutine load_ppn
         integer :: i, fh, cnt
         character(len=5) :: iso
         real(r8) :: x
         open( file = ini_filename, newunit = fh )
         do i = 1, 7
            read(fh,*) ! skip header
         end do
         cnt = 1
         do while ( .true. )
            read( fh, ppn_fmt, end = 55 ) x, iso
            if ( x > 1.e-60_r8 ) then
               iso_in(cnt) = iso; x_in(cnt) = x
               cnt = cnt + 1
            end if
            if (cnt > nsp) stop 'load_abundances: the network is missing isotopes for which you want to initialise the abundances'
         enddo
         55 num_iso_in = cnt - 1
         close( fh )
   end subroutine load_ppn

   !---
   !> @brief read abundances from a frischknecht-format ascii file
   !---
   subroutine load_urs_frischknecht_xin
         integer :: fh, cnt
         character(len=5) :: iso
         real(r8) :: x
         fh = get_unused_file_handle()
         cnt = 1
         open( fh, file = ini_filename )
         do while ( .true. )
            read( fh, frisch_fmt, end = 55 ) iso, x
            if ( x > 1.e-60_r8 ) then
               iso_in(cnt) = iso; x_in(cnt) = x
               cnt = cnt + 1
            end if
            if (cnt > nsp) stop 'load_abundances: the network is missing isotopes for which you want to initialise the abundances'
         enddo
         55 num_iso_in = cnt - 1
         close( fh )
   end subroutine load_urs_frischknecht_xin

   !---
   !> @brief read abundances from a ascii file written with abu_from_masscoordinate
   !---
   subroutine load_abu_from_masscoordinate
         integer :: fh, cnt, iso_a
         character(len=5) :: cdummy, iso
         character(len=2) :: iso_elname
         real(r8) :: x
         fh = get_unused_file_handle()
         open ( fh, file = ini_filename )
         cnt = 1
         do while ( .true. )
            read( fh, '(a1)', end=9344) cdummy
            if ( cdummy == 'D' ) then
               backspace( fh )
               read( fh, * ) cdummy, iso_elname, iso_a, x
               if ( x > 1.e-60_r8 ) then
                  write( iso_in(cnt), '(A2,I3)' ) iso_elname, iso_a
                  x_in(cnt) = x
                  cnt       = cnt + 1
               end if
            endif
         end do
         close( fh )
         9344 num_iso_in = cnt - 1
   end subroutine load_abu_from_masscoordinate

   !---
   !> @brief read abundances from a ascii file written with nugridpy
   !---
   subroutine load_nugridpy_utils_xin
         integer :: fh, cnt
         character(len=80) :: cheader
         fh = get_unused_file_handle()
         open ( fh, file = ini_filename )
         read ( fh, '(1x,A)') cheader
         write (*,*) 'Reading initial abundance with the following header'
         write (*,*) 'comment:'
         write (*,*) cheader
         read( fh, '(1x,A)') cheader

         cnt = 1
         do while ( .true. )
            read( fh, nugrid_fmt, end = 56) iso_in(cnt), x_in(cnt)
            cnt = cnt + 1
            if (cnt > nsp) stop 'load_abundances: the network is missing isotopes for which you want to initialise the abundances'
         enddo
         56 num_iso_in = cnt - 1
         close( fh )
   end subroutine load_nugridpy_utils_xin

   !---
   !> @brief read abundances from a plainascii file with isotope name, abundance
   !---
   subroutine load_plain_ascii
         integer :: fh, cnt
         fh = get_unused_file_handle()
         open ( fh, file = ini_filename )
         write (*,*) 'Reading initial abundances from ', ini_filename

         cnt = 1
         do while ( .true. )
            read( fh, "(a5,es23.15)", end = 56) iso_in(cnt), x_in(cnt)
            cnt = cnt + 1
            if (cnt > nsp) stop 'load_abundances: the network is missing isotopes for which you want to initialise the abundances'
         enddo
         56 num_iso_in = cnt - 1
         close( fh )
   end subroutine load_plain_ascii


end module abundances
