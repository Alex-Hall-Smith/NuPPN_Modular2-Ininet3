!> @author Marco Pignatari, Falk Heriwg, Samuel Jones
!
!> @brief This module handles the reading and writing of the networksetup.txt files, which
!> contain the details of which species and reactions are to be included, and from where the
!> reaction rates are to be sourced
!
!> @details Detailed description of the purpose of this module to appear here...
!
module networksetup
   use nuc_data, only: zis, considerisotope, considerreaction, niso, isomeric_state
   use rates, only: k1, k2, k3, k4, k5, k6, k7, k8, v, ant, znt, znf, rates_hash_locations_for_merge
   use utils, only: r8, i4
   use array_sizes
   use constants, only: ZERO
   use reaction_info, only: lab, labb, ilabb, bind_energy_diff, rfac
   use communication

   implicit none
   private

   ! file names and handles
   integer, parameter :: nwsetup_fh = 77
   character(len=16), parameter :: nwsetup_fn  = "networksetup.txt"
   character(len=17), parameter :: nwsetup2_fn = "networksetup2.txt"

   ! formatting
   character(len=256), parameter :: &
         fmtstars          = "(72('*'))", &
         fmtspeciesheader1 = "(8X,'TABLE OF ISOTOPES ',e10.3,1x,e10.3)", &
         fmtspeciesheader2 = "(2X,4X,'  NO.',1X,'SYMBOL T/F MASSNR PROTON')", &
         fmtspecies        = "(7X,I5,3X,A5,1X,L5,2X,0PF4.0,3X,F4.0,3x,I1)", &
         fmtreactionheader = "(10X,'REACTION NETWORK',35X,'V(i) Nasv(*rho)')", &
         fmtreaction       = "(2X,I5,1X,L1,1X,I2,2X,A5,'  + ',I2,2X,A5,2X,'->',1X,I2,2X,A5,'  + ',"//&
                             "I2,2X,A5,3X,ES9.3,2x,A5,2x,A5,2x,I2,2x,ES10.3,2x,ES10.3)", &
         fmtreactionread   = "(2X,I5,1X,L1,1X,I2,2X,A5,4X,I2,2X,A5,2X,2X,1X,I2,2X,A5,4X,"//&
                             "I2,2X,A5,3X,ES9.3,2x,A5,2x,A5,2x,I2,2x,ES10.3,2x,ES10.3)"

   ! from common/nwsetup1
   integer(i4), public :: &
         inumer, &  !< number of isomers in networksetup.txt file
         inumerr, & !< number of reactions involviing isomers in networksetup.txt file
         nvnc1,  &  !< number of isotopes in networksetup.txt file, including isomers
         nrnc1      !< number of reactions in networksetup.txt file, including ones with isomers

   public read_networksetup, write_networksetup, write_networksetup2
#ifndef PPN
   public networksetup_broadcasts
#endif

contains

   !---
   !> @brief Write the file networksetup.txt
   !---
   subroutine write_networksetup(t9, rho, ngis, nrcp, ngir, an, zn)
         integer :: i, ngis, nrcp, ngir
         real(r8) :: t9, rho, an(nsp), zn(nsp)

         open (nwsetup_fh, file = nwsetup_fn)
         call output_network(t9, rho, ngis, nrcp, ngir, an, zn)
         close(nwsetup_fh)
   end subroutine write_networksetup

   !---
   !> @brief Write the file networksetup2.txt
   !---
   subroutine write_networksetup2(t9, rho, ngis, nrcp, ngir, an, zn)
         integer :: i, ngis, nrcp, ngir
         real(r8) :: t9, rho, an(nsp), zn(nsp)

         open (nwsetup_fh, file = nwsetup2_fn)
         call output_network(t9, rho, ngis, nrcp, ngir, an, zn)
         close(nwsetup_fh)
   end subroutine write_networksetup2


   !> output network details: isotopes, reactions etc to a file
   subroutine output_network(t9, rho, ngis, nrcp, ngir, an, zn)
         use nuc_data, only: isomeric_state
         integer :: i, &
               ngis, & !< number of species including isomers ?
               nrcp, & !< number of  "charged particle" reactions (i.e. reactions from VITAL)
               ngir    !< total number of reactions including isomers (I think)
         real(r8) :: &
               t9, &      !< temperature (GK)
               rho, &     !< density (gcc)
               an(nsp), & !< mass numbers @todo can't this be taken from nuc_data?
               zn(nsp)    !< proton numbers @todo can't this be taken from nuc_data?

         write(nwsetup_fh, fmtstars         )
         write(nwsetup_fh, fmtspeciesheader1) t9, rho
         write(nwsetup_fh, *                )
         write(nwsetup_fh, fmtspeciesheader2)
         write(nwsetup_fh, *                ) ngis,' =  NGIS', inumer,' =  ISOMERS'
         
         ! species
         !--------
         do i = 1, ngis
            write(nwsetup_fh, fmtspecies) i, zis(i), considerisotope(i), an(i), zn(i), isomeric_state(i)
         end do

         write(nwsetup_fh, fmtstars         )
         write(nwsetup_fh, *                )
         write(nwsetup_fh, fmtreactionheader)
         write(nwsetup_fh, *                )

         ! "charged particle" reactions (i.e. from VITAL I guess)
         !-------------------------------------------------------
         write(nwsetup_fh, *                ) ngir, ' =  NGIR', inumerr, ' = + REAZ ISOMERS'

         do i = 1, nrcp
            write(nwsetup_fh, fmtreaction   ) i, considerreaction(i), &
                  k2(i), zis(k1(i)), k4(i), zis(k3(i)), k8(i), zis(k7(i)), k6(i), zis(k5(i)), &
                  v(i), lab(i), labb(i), ilabb(i), rfac(i), bind_energy_diff(i)
         end do

         ! "other", sometimes called "neutron capture" reactions
         !------------------------------------------------------
         ! (but I think it just means non-VITAL)
         do i = nrcp + 1, ngir
            write(nwsetup_fh, fmtreaction   ) i, considerreaction(i), &
                  k2(i), zis(k1(i)), k4(i), zis(k3(i)), k8(i), zis(k7(i)), k6(i), zis(k5(i)), &
                  v(i), lab(i), labb(i), ilabb(i), rfac(i), bind_energy_diff(i)
         end do
   end subroutine output_network


   !!---
   !!> @brief Writes file networksetup2.txt
   !!---
   !subroutine write_networksetup2(nvnc, nrnc, an, zn)
         !integer :: i, &
               !nvnc, & !< number of 'charged particle' reactions (i.e. from VITAL)
               !nrnc    !< number of 'neutron capture' reactions (i.e. non-VITAL)
         !real(r8) :: &
               !an(nsp), & !< mass numbers
               !zn(nsp)    !< proton numbers

         !open(nwsetup_fh, file = nwsetup2_fn)

         !write(nwsetup_fh, fmtstars          )
         !write(nwsetup_fh, *                 ) "        TABLE OF ISOTOPES"
         !write(nwsetup_fh, fmtspeciesheader2 )
         !write(nwsetup_fh, *                 ) nvnc, ' =  NVNC'

         !do i = 1, nvnc
            !write(nwsetup_fh, fmtspecies     ) i, zis(i), considerisotope(i), an(i), zn(i), isomeric_state(i)
         !end do

         !write(nwsetup_fh, fmtstars          )
         !write(nwsetup_fh, fmtreactionheader )
         !write(nwsetup_fh, *                 ) nrnc, ' =  NRNC'

         !do i = 1, nrnc
            !write(nwsetup_fh, fmtreaction    ) i, considerreaction(i), &
                  !k2(i), zis(k1(i)), k4(i), zis(k3(i)), k8(i), zis(k7(i)), k6(i), zis(k5(i)), &
                  !v(i),  lab(i), labb(i), ilabb(i), rfac(i)
         !end do

         !close(nwsetup_fh)
   !end subroutine write_networksetup2


   !---
   !> @brief Read file networksetup.txt and initialise network that way
   !---
   subroutine read_networksetup(an, zn, t9_nw_ini, rho_nw_ini)
         use nuc_data, only: ispe
         use communication
         integer :: i, l, nnnr, nnnz
         real(r8) :: an(nsp), zn(nsp), t9_nw_ini, rho_nw_ini
         character(len=5) :: spe1, spe2, spe3, spe4
         character(len=256) :: tmp1
         character(len=25) :: tmp2

         bind_energy_diff(:) = ZERO

         if (master) then
            open(nwsetup_fh, file = nwsetup_fn)

            read(nwsetup_fh, *               )
            read(nwsetup_fh, "(a25,2(e11.3))") tmp2, t9_nw_ini, rho_nw_ini
            read(nwsetup_fh, *               )
            read(nwsetup_fh, *               )
            read(nwsetup_fh, '(7x,I5,16x,I5)') nvnc1,inumer

            do i = 1, nvnc1
               read ( nwsetup_fh, fmtspecies ) &
                     nnnz, zis(i), considerisotope(i), an(i), zn(i), isomeric_state(i)
               ! ^_^ failsafe for historic bug with fake isotopes
               if ( zis(i) /= 'NEUT ' .and. ( an(i) == 1 .and. zn(i) == 0 ) ) then
                  stop "Duplicate neutrons in networksetup.txt"
               end if
            end do

            ! *** creation of index array niso
            ! *** third argument is 1 (ground) o 2 (isomeric).
            ! *** It will account of the ground/isomeric state.
            niso(:,:,:) = 0
            do i = 1, nvnc1
               niso(nint(an(i)),nint(zn(i)),isomeric_state(i)) = i
            end do

            read(nwsetup_fh, *               )
            read(nwsetup_fh, *               )
            read(nwsetup_fh, *               )
            read(nwsetup_fh, *               )
            read(nwsetup_fh, '(7x,I5,16x,I5)') nrnc1, inumerr

            do i = 1, nrnc1
               read( nwsetup_fh, fmtreactionread  ) nnnr, considerreaction(i), &
                     k2(i), spe1, k4(i), spe2, k8(i), spe3, k6(i), spe4, &
                     v(i), lab(i), labb(i), ilabb(i), rfac(i), bind_energy_diff(i)

               k1(i) = ispe(spe1)
               k3(i) = ispe(spe2)
               k5(i) = ispe(spe4)
               k7(i) = ispe(spe3)
            end do

            close(nwsetup_fh)

            ! make sure that all reactions from networksetup.txt including species
            ! that are set as false by hand in networksetup.txt are set to .false.
            do i=1,nvnc1
               if ( .not. considerisotope(i) .and. ispe(zis(i)) .ne. ispe("OOOOO") ) then
                  write(*,*) 'rnetw2008: false species=',zis(i)
                  do l = 1, nrnc1
                     if ( ispe(zis(i)) == k1(l) .or. ispe(zis(i)) == k3(l) ) then
                        considerreaction(l) = .false.
                     end if
                     if ( ispe(zis(i)) == k5(l) .or. ispe(zis(i)) == k7(l) ) then
                        considerreaction(l) = .false.
                     end if
                  end do
               end if
            end do
         end if ! master

#ifndef PPN
         call networksetup_broadcasts()
#endif

! creation of ant,znt. Used in rnetw2008 to merge rates in ppn network.
         do i = 1, nvnc1
            do l = 1, nrnc1
               if ( zis(i) == zis(k1(l)) ) then
                  ant(l) = int(an(i))
                  znt(l) = int(zn(i))
               end if
               if ( zis(i) == zis(k3(l)) )then
                  znf(l) = int(zn(i))
               end if
            end do
         end do

         ! hash the locations of each rate source in the main rate array
         call rates_hash_locations_for_merge(nrnc1, ilabb)

   end subroutine read_networksetup


#ifndef PPN
   subroutine networksetup_broadcasts
      real(r8) :: an(nsp), zn(nsp)
      common / cnetw    / an, zn
      call broadcast(k1); call broadcast(k2); call broadcast(k3); call broadcast(k4)
      call broadcast(k5); call broadcast(k6); call broadcast(k7); call broadcast(k8)
      call broadcast(bind_energy_diff); call broadcast(rfac); call broadcast(ilabb)
      call broadcast_ch_arr(lab); call broadcast(v); call broadcast(considerreaction)
      call broadcast_ch_arr(labb)
      call broadcast(nrnc1); call broadcast(inumerr)
      call broadcast(niso)
      call broadcast(nvnc1); call broadcast(inumer)
      call broadcast_ch_arr(zis); call broadcast(considerisotope); call broadcast(an)
      call broadcast(zn); call broadcast(isomeric_state)
   end subroutine networksetup_broadcasts
#endif

end module networksetup
