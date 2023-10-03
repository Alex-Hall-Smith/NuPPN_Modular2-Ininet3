! MODULE: NKK04
!
!> @author Samuel Jones
!
!> @brief This module provides weak reaction rates from J-U. Nabi, Klapdor-Kleingrothaus, and H. Volker. ADNDT, 88 2 (2004) 237-476.
!
!> @details NKK04 rate coverage overlaps with several other rate sources,

!> so there are some options for how to merge the rates:
!> nkk_mode = 1 : don't use NKK rates at all
!> nkk_mode = 2 : use NKK rates everywhere except where LMP00 rates are available
!> nkk_mode = 3 : use NKK rates everywhere they are available
!
module nkk04
   use frame_knobs, only: iolevel
   use utils, only: r8, get_unused_file_handle, irate, idrdt, idrdd, generic_interp_2d
   use physics_knobs, only: nkk_mode
   use constants   
#ifndef PPN
   use communication
#endif
   implicit none
   private

   character(len=256), parameter :: nkk_fn = '../NPDATA/nkk04_mmc1_fixed.txt'
   integer, parameter :: &
         num_nkk_rates = 1304, &
         nT = 12, & !< temperature (T9) dimension
         nr = 11, & !< electron density dimension
         nkkminA = 40, nkkmaxA = 100, nkkminZ = 19, nkkmaxZ = 51
   character(len=4) :: weaktype(num_nkk_rates) ! 'beta' or 'ecap'
   real(r8) :: &
         T9(nT), &                           !< temperature coordinates
         log10T(nT), &                       !< log10 temperature coordinates
         lyr(nr), &                          !< electron density coordinates
         raw_rates(num_nkk_rates,nr,nT), &   !< the whole grid of reaction rates
         nkk_rate(num_nkk_rates), &          !< holds the rate for each reaction, interpolated to the required T and rho
         dlT, dT9, dlyr, &                   !< grid spacings
         minT9, maxT9, minlyr, maxlyr, minlT, maxlT   !< grid bounds
   integer :: &
         Ai(num_nkk_rates), & !< parent nucleus A for each reaction
         Zi(num_nkk_rates), & !< parent nucleus Z for each reaction
         Af(num_nkk_rates), & !< daughter nucleus A for each reaction
         Zf(num_nkk_rates), & !< daughter nucleus Z for each reaction
         nkk_ntrans( nkkminA:nkkmaxA, nkkminZ:nkkmaxZ, 1, 13:14 ) !< (A, Z, 'state', reaction type); for merging into main network
                                                      ! reaction type 13 = b- or e+-capture; 14 = b+ or e--capture

   public nkk_init, nkk_rate, nkkminA, nkkmaxA, nkkminZ, nkkmaxZ, nkk_ntrans, nkk_interpolate_rates

contains

   !---
   !> @brief Initialise NKK04 module
   !---
   subroutine nkk_init()
         integer :: i

         if (nkk_mode == 1) return

         nkk_ntrans(:,:,:,:) = 0
#ifndef PPN
         if (ipid == 0) then
#endif
            call nkk_read_data()
#ifndef PPN
         end if
         call nkk_broadcasts
#endif
         do i = 1, num_nkk_rates
            if ( Zf(i) > Zi(i) ) then
               weaktype(i) = 'beta'
               nkk_ntrans( Ai(i), Zi(i), 1, 13 ) = i
            else
               weaktype(i) = 'ecap'
               nkk_ntrans( Ai(i), Zi(i), 1, 14 ) = i
            end if
         end do

         dlyr = 1._r8
         minT9    = minval(T9);  maxT9    = maxval(T9)
         minlyr   = minval(lyr); maxlyr   = maxval(lyr)

   end subroutine nkk_init


   !---
   !> @brief Perform bilinear interpolation (in log10 electron density and temperature in GK)
   !> of NKKreaction rates
   !---
   subroutine nkk_interpolate_rates(rho_in, ye_in, T9_in)
         real(r8), intent(in) :: &
               rho_in, & !< density at which the rates are needed
               ye_in, &  !< electron fraction at which the rates are needed
               T9_in     !< temperature (GK) at which the rates are needed
         real(r8) :: lyr_in, myT9, mylyr, &
               xa, xb, ya, yb, &       !< interpolation variables
               T91, T92, lyr1, lyr2, & !< bounding values of the interpolation
               f11, f21, f12, f22, &   !< bounding data values
               fA, fB, fC, fD
         real(r8), dimension(num_nkk_rates) :: lr, dlrdt9, dlrdlyr
         integer :: i, &
               idxlT, idxT9, idxlyr           !< indices of lower bounding values

         lyr_in = log10( rho_in * ye_in )
         ! clip values to table edge
         myT9  = min( maxT9,  max( T9_in,  minT9  ) )
         mylyr = min( maxlyr, max( lyr_in, minlyr ) )

         ! ^_^ get index using grid spacing (constant spacing in log(rho * Ye) )
         idxT9    = minloc( abs( T9(:) - myT9 ), dim = 1 )
         if ( T9(idxT9) > myT9 ) idxT9 = idxT9 - 1
         idxlyr   = floor( ( mylyr - minlyr ) / dlyr ) + 1

         ! interpolants (constant across rates)
         T91   = T9  ( idxT9    )
         T92   = T9  ( idxT9+1  )
         lyr1  = lyr ( idxlyr   )
         lyr2  = lyr ( idxlyr+1 )

         dT9 = T92 - T91

         xb = ( mylyr - lyr1 ) / dlyr
         xa = 1._r8 - xb
         yb = ( myT9 - T91 ) / dT9
         ya = 1._r8 - yb

         call generic_interp_2d(nr, nT, lyr, T9, num_nkk_rates, raw_rates(:,:,:), mylyr, myT9, &
               idxlyr, idxlyr+1, idxT9, idxT9+1, lr, dlrdt9, dlrdlyr)

         nkk_rate(:) = 10._r8 ** lr(:)

         if (.false.) then
            print *, '************************************'
            print *, 'you asked for rate at T9, logYeRho =',T9_in, lyr_in, myT9, mylyr
            print *, 'this was found to be between temperatures', T9(idxT9), T9(idxT9+1) 
            print *, 'and densities', lyr(idxlyr), lyr(idxlyr+1)
            print *, 'with indices', idxT9, idxT9+1, idxlyr, idxlyr+1
            print *, 'bounding values of last rate:', f11, f21, f12, f22
            print *, 'interpolation constants:', xa, xb, ya, yb, fA, fB, fC, fD
            print *, 'table limits:', minT9, maxT9, minlyr, maxlyr
            print *, '************************************'
         end if

   end subroutine nkk_interpolate_rates


   !---
   !> @brief read the NKK reaction rate data tables
   !---
   subroutine nkk_read_data()
         integer :: nkk_fh
         integer :: i, j, k, l, line, ierr, Tdim, rdim
         real(r8) :: rate1, rate2
         ! ^_^ scratch variables for dummy reading
         integer :: idum
         character(len=2) :: cdum
         real(r8) :: rdum, rdum_arr(4)

         nkk_fh = get_unused_file_handle()
         if ( iolevel >= 1) print *, 'reading weak rates from NKK04...'
         open( nkk_fh, file = nkk_fn, status='old', action = 'read' )

         do i = 1, num_nkk_rates
            ! ^_^ read a rate
            read (nkk_fh, *) Ai(i), idum
            read (nkk_fh, *) cdum, cdum
            read (nkk_fh, *) Zi(i), idum, Zf(i), idum
            ! ^_^ skip 5 lines (blanks plus column headers)
            do j = 1, 5
               read (nkk_fh, *)
            end do
            ! ^_^ read coordinates and rate data (sum primary and complimentary rates)
            do j = 1, nr
               do k = 1, nT
                  read (nkk_fh,*) lyr(j), T9(k), rdum, rate1, rate2, (rdum_arr(l), l = 1, 4)
                  raw_rates(i,j,k) = log10(10._r8 ** rate1 + 10._r8 ** rate2)
               end do
            end do
            read (nkk_fh, *) ! blank line
         end do

         close( nkk_fh )
   end subroutine nkk_read_data



#ifndef PPN
   !---
   !> @brief Broadcast NKK rate data from master to slaves in parallel frames
   !---
   subroutine nkk_broadcasts
         call broadcast(raw_rates)
         call broadcast(Zi)
         call broadcast(Zf)
         call broadcast(Ai)
         call broadcast(Af)
         call broadcast(lyr)
         call broadcast(T9)
         return
   end subroutine nkk_broadcasts
#endif


end module nkk04
