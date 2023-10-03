! MODULE: alpha_decays
!
!> @author Marco Pignatari, Samuel Jones
!
!> @brief This module provides alpha-decay half lives from the Karlsruhe Nuclidkarte
!
!
module alpha_decays
   use constants
   use nuc_data, only: ispe
   use utils, only: r8, irate, idrdt, idrdd, get_unused_file_handle
   use communication
   implicit none
   private

   character(len=256), parameter :: adcy_fn = '../NPDATA/karlsruhe_nuclidkarte_2006.dat'
   character(len=5), parameter :: adcy_name = 'MPG06' !< tag indicating source of rates
   integer, parameter :: num_adcy = 331 !< number of species for which there are measurements
   integer :: &
         num_packed_adcy, & !< number of alpha decays that will actually be used in the network
         iadcy, adcy_fh
   logical :: &
         alpha_decay_mask(num_adcy), & !< true where parent and daughter are in network
         adcy_masks_exist              !< whether or not alpha_decay_mask array has been created
   real(r8) :: &
         adcy_half_life(num_adcy), &   !< alpha decay half lives
         adcy_rate(num_adcy), &        !< alpha decay rates \f$ (\frac{\ln 2}{\tau_{1/2}}) \f$
         adcy_rate_packed(num_adcy)    !< alpha decay rates \f$ (\frac{\ln 2}{\tau_{1/2}}) \f$, packed
   character(len=4) :: &
         units(num_adcy)            !< time unit for the decay half life
   character(len=5) :: &
         p_iso(num_adcy), &         !< parent isotope
         d_iso(num_adcy), &         !< daughter isotope
         p_iso_packed(num_adcy), &  !< parent isotope (packed)
         d_iso_packed(num_adcy)     !< daughter isotope (packed)

   public adcy_name, alpha_decays_init, adcy_masks_exist, adcy_rate_packed, p_iso_packed, d_iso_packed, &
#ifndef PPN
         alpha_decays_broadcasts, &
#endif
         num_packed_adcy, iadcy, alpha_decays_create_masks


contains

   !---
   !> @brief Initialises alpha decays module
   !---
   subroutine alpha_decays_init()

         adcy_masks_exist = .false.
         iadcy = -1

         if (master) then
            call read_alpha_decay_data()
         end if

#ifndef PPN
         call broadcast_ch_arr(p_iso)
         call broadcast_ch_arr(d_iso)
         call broadcast(adcy_half_life)
         call broadcast(adcy_rate)
#endif

   end subroutine alpha_decays_init


   !---
   !> @brief read alpha decay half lives
   !---
   subroutine read_alpha_decay_data()

         integer :: z, a, i

         adcy_fh = get_unused_file_handle()

         open( unit = adcy_fh, file = adcy_fn, status = 'old' )

         adcy_half_life(:) = 0._r8
         adcy_rate(:) = 0._r8

         ! header
         read(adcy_fh,*)
         read(adcy_fh,*)       

         ! data
         do i = 1, num_adcy
            read(adcy_fh,'(2x,a5,4x,a5,2x,i3,6x,i3,3x,e9.2,4x,a4)') &
                  p_iso(i), d_iso(i), z, a, adcy_half_life(i), units(i)
            select case(units(i))
            case('muse')
               adcy_half_life(i) = adcy_half_life(i) * mus2sec
            case('msec')
               adcy_half_life(i) = adcy_half_life(i) * msec2sec
            case('secs')
               continue
            case('mins')
               adcy_half_life(i) = adcy_half_life(i) * min2sec
            case('hour')
               adcy_half_life(i) = adcy_half_life(i) * hrs2sec
            case('days')
               adcy_half_life(i) = adcy_half_life(i) * day2sec
            case('year')
               adcy_half_life(i) = adcy_half_life(i) * yrs2sec
            end select
            adcy_rate(i) = ln2 / adcy_half_life(i)
         end do

         close( adcy_fh )

   end subroutine read_alpha_decay_data


   !---
   !> @brief check if both parent and daughter of alpha decay are included in current network,
   !> and create an array of just these parent and daughter nuclei and the decay rate
   !---
   subroutine alpha_decays_create_masks()
         integer :: i

         adcy_rate_packed  (:) = 0._r8
         p_iso_packed      (:) = '     '
         d_iso_packed      (:) = '     '

         num_packed_adcy = 0

         do i = 1, num_adcy
            if (ispe(p_iso(i)) .ne. -1 .and. ispe(d_iso(i)) .ne. -1) then
               num_packed_adcy      = num_packed_adcy + 1
               alpha_decay_mask(i)  = .true.
               adcy_rate_packed (num_packed_adcy) = adcy_rate(i)
               p_iso_packed     (num_packed_adcy) = p_iso(i)
               d_iso_packed     (num_packed_adcy) = d_iso(i)
            end if
         end do

         adcy_masks_exist = .true.

   end subroutine alpha_decays_create_masks

#ifndef PPN
   !---
   !> @brief Broadcast alpha decay data from master to slaves in parallel frames
   !---
   subroutine alpha_decays_broadcasts()
         call broadcast(alpha_decay_mask)
         call broadcast(num_packed_adcy)
         call broadcast(adcy_rate_packed)
         call broadcast_ch_arr(p_iso_packed)
         call broadcast_ch_arr(d_iso_packed)
         call broadcast(iadcy)
   end subroutine alpha_decays_broadcasts
#endif

end module alpha_decays
