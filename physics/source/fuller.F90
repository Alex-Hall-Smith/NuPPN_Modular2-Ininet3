module fuller
   use utils, only : r8, get_unused_file_handle, irate, idrdt, idrdd, generic_interp_2d
   use array_sizes
   use constants
#ifndef PPN
   use communication
#endif

   implicit none

   private

   ! ^_^ file name
   character(len=256), parameter :: fuller_fn = '../NPDATA/fullerweak.txt'

   ! ^_^ table dimensions and independent variables (coordinates)
   integer, parameter :: &
         fuller_vdim = 189 * 2, &
         fuller_densdim = 11, &
         fuller_t9dim = 13

   real(r8) :: fuller_t9coo(fuller_t9dim), fuller_rhocoo(fuller_densdim)
   ! Z and A of "neg. daughter" species
   integer, dimension(fuller_vdim)  :: fuller_mz, fuller_ma
   ! Q-value
   real(r8), dimension(fuller_vdim) :: fuller_dmf
   ! rate requested T9 and log10(rho*ye)
   real(r8) :: fuller_rates(fuller_vdim)
   ! ^_^ for merging in ppn
   integer  :: fuller_ntrans(60,0:30,1,13:14)
   ! ^_^ per-reaction tabulated data
   real(r8), dimension(fuller_t9dim,fuller_densdim) :: ff1   ! fermi energy
   real(r8), dimension(fuller_vdim,fuller_t9dim,fuller_densdim) :: &
         ff2,  & ! lbeta(+/-)
         ff3,  & ! leps(-/+)
         ff4,  & ! lsum(+/-)
         ff5     ! lnu/lnubar

   ! ^_^ publish
   public fuller_init, fuller_interpolate_rates, fuller_ntrans, fuller_vdim, fuller_rates

contains

   subroutine fuller_init()
         implicit none

         fuller_ntrans(:,:,:,:) = 0

         ! definition of temperature and density coordinates in fuller tables
         fuller_t9coo(1)    =    log10(0.01_r8)
         fuller_t9coo(2)    =    log10(0.10_r8)
         fuller_t9coo(3)    =    log10(0.20_r8)
         fuller_t9coo(4)    =    log10(0.40_r8)
         fuller_t9coo(5)    =    log10(0.70_r8)
         fuller_t9coo(6)    =    log10(1.00_r8)
         fuller_t9coo(7)    =    log10(1.50_r8)
         fuller_t9coo(8)    =    log10(2.00_r8)
         fuller_t9coo(9)    =    log10(3.00_r8)
         fuller_t9coo(10)   =    log10(5.00_r8)
         fuller_t9coo(11)   =   log10(10.00_r8)
         fuller_t9coo(12)   =   log10(30.00_r8)
         fuller_t9coo(13)   =  log10(100.00_r8)

         fuller_rhocoo(1)   =  1._r8
         fuller_rhocoo(2)   =  2._r8
         fuller_rhocoo(3)   =  3._r8
         fuller_rhocoo(4)   =  4._r8
         fuller_rhocoo(5)   =  5._r8
         fuller_rhocoo(6)   =  6._r8
         fuller_rhocoo(7)   =  7._r8
         fuller_rhocoo(8)   =  8._r8
         fuller_rhocoo(9)   =  9._r8
         fuller_rhocoo(10)  = 10._r8
         fuller_rhocoo(11)  = 11._r8


#ifndef PPN
         if ( ipid == 0 ) then

            call fuller_read_data()

         end if

         call fuller_broadcasts()
#else
         call fuller_read_data()
#endif
   end subroutine fuller_init



   subroutine fuller_interpolate_rates( t9in, rhoin, yein )
         use utils, only: bsearch_r8, r8
         implicit none
         real(r8), intent(in) :: t9in, rhoin, yein
         integer :: i, it90_1, it90_2, irho_1, irho_2
         real(r8) :: t90, t90_1, t90_2, rho0, rho_1, rho_2, xa, xb, ya, yb
         real(r8), dimension(fuller_vdim) :: ff11, ff12, ff21, ff22, lr, dlrdt9, dlrdlyr

         ! use log10(t9) instead of t9 in the interpolations
         t90  = log10(t9in)
         rho0 = log10(rhoin * yein)

         ! clip the temperature and density to the table
         t90   = min(max(t90, fuller_t9coo(1)), fuller_t9coo(fuller_t9dim))
         rho0  = min(max(rho0, fuller_rhocoo(1)), fuller_rhocoo(fuller_densdim))

         ! bisection search
         call bsearch_r8( fuller_t9coo, t90, it90_1 )
         call bsearch_r8( fuller_rhocoo, rho0, irho_1 )
         it90_2 = it90_1 + 1; irho_2 = irho_1 + 1
         t90_1 = fuller_t9coo(it90_1); t90_2 = fuller_t9coo(it90_2)
         rho_1 = fuller_rhocoo(irho_1); rho_2 = fuller_rhocoo(irho_2)

         ! ^_^ bilinear interpolation of ff4(:)

         call generic_interp_2d(fuller_t9dim, fuller_densdim, fuller_t9coo, fuller_rhocoo, fuller_vdim, &
               ff4(:,:,:), t90, rho0, it90_1, it90_2, irho_1, irho_2, lr, dlrdt9, dlrdlyr)

         fuller_rates(:)    = 10._r8 ** lr(:)

   end subroutine fuller_interpolate_rates



   subroutine fuller_read_data()

         ! *** note: the third argument of fuller_ntrans  is 1, because in fuller
         !           table there is just thermalized species (possible problem,
         !           e.g.Al26).

         integer     :: i, j, jj, nz1, na1, na2, nz2, fuller_fh
         real(r8)    :: dm1, dm2
         character   :: cdummy

         fuller_fh = get_unused_file_handle()
         open(unit = fuller_fh, file = fuller_fn, status = 'old')

         i = 0

         27 read(fuller_fh, '(6x,a1)', end = 992) cdummy

         if (cdummy == '#') goto 992

         backspace(fuller_fh)

         read(fuller_fh, '(23x,i2,11x,i3,5x,f8.4)') nz1, na1, dm1

         fuller_mz(i+1)  = nz1
         fuller_ma(i+1)  = na1
         fuller_dmf(i+1) = dm1
         fuller_ntrans(na1,nz1,1,14) = i + 1

         read(fuller_fh,'(23x,i2,11x,i3,5x,f8.4)') nz2, na2, dm2

         fuller_mz(i+2)  = nz2
         fuller_ma(i+2)  = na2
         fuller_dmf(i+2) = dm2
         fuller_ntrans(na2,nz2,1,13) = i + 2

         read(fuller_fh,*)
         read(fuller_fh,*)

         ! ^_^ rates

         do j = 1, fuller_densdim
            do jj = 1, fuller_t9dim
               read(fuller_fh,'(10x,2(f7.3),7(1x,f7.3))') ff1(jj,j), &
                     ff2(i+1,jj,j),ff3(i+1,jj,j),ff4(i+1,jj,j),ff5(i+1,jj,j), &
                     ff2(i+2,jj,j),ff3(i+2,jj,j),ff4(i+2,jj,j),ff5(i+2,jj,j)
            end do
         end do

         read(fuller_fh,*)

         i = i + 2

         goto 27

         992 continue

         close(fuller_fh)

         return

   end subroutine fuller_read_data


#ifndef PPN
   subroutine fuller_broadcasts()

         call broadcast(fuller_mz)
         call broadcast(fuller_ma)
         call broadcast(fuller_dmf)
         call broadcast(fuller_ntrans)
         call broadcast(ff1)
         call broadcast(ff2)
         call broadcast(ff3)
         call broadcast(ff4)
         call broadcast(ff5)

         return

   end subroutine fuller_broadcasts
#endif




end module fuller
