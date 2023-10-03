module nuc_data
   use array_sizes, only: nsp, nre
   use utils, only: int_hash_table, char_hash_table, create_hash_table, hash, Cantor_hash, r8, i4
   implicit none
   private

   character(len=5) :: zis(nsp)           !< isotope names
   logical :: considerisotope(nsp)        !< switches for adaptive network
   logical :: considerreaction(nre)       !< switches for adaptive network

   ! isotope database
   character(len=256), parameter :: isotopedatabase_fn = 'isotopedatabase.txt'
   integer, parameter :: isotopedatabase_fh = 3
   integer, dimension(nsp) :: atot, ztot, ntot
   real(r8) :: amin(nsp) ! why??
   character(len=2) :: eltot(nsp)

   ! destroy when possible:
   common/networklist1/ztot
   common/networklist2/amin,atot,eltot

   real(r8), dimension(nsp) :: &
         anetw, &  !< mass number
         znetw, &  !< proton number
         t_half    !< half lives
   common/cnetw/anetw,znetw

   character :: t_half_unit(nsp)


   !> names of elements:
   character(len=2) :: element_name_lower(0:98), element_name_upper(0:98)

   !> pointer to isotope in global state arrays given A, Z and isomeric state
   integer(i4), allocatable, public :: niso(:,:,:)

   !> state of each species (1, 2, 3, ...), 1 is ground (thermalised?), 2 is isomer
   integer(i4), public :: isomeric_state(nsp)

   !> for hashing
   type(int_hash_table)    :: zis_name_to_idx
   type(char_hash_table)   :: az_to_zis_hash_idx

   public zis, considerisotope, considerreaction, zis_name_to_idx, &
         zis_create_hash, create_az_to_zis_hash, nuc_data_init, ispe, epsi, t_half, t_half_unit, &
         atot, ztot, ntot, eltot, amin, read_isotopedatabase, atomic

contains

   subroutine read_isotopedatabase(isotot, considerisoinlist)
         use communication
         use utils
         integer :: i, isotot
         logical :: considerisoinlist(nsp)

         if (master) then
            print *, "init nuc data ... "
            open(unit = isotopedatabase_fh, file = isotopedatabase_fn, status = 'old')

            read(isotopedatabase_fh,*)
            read(isotopedatabase_fh,*)
            read(isotopedatabase_fh,*)

            do i = 1, isotot
               read(isotopedatabase_fh, *, end=98778) ztot(i), atot(i), eltot(i), amin(i), &
                     considerisoinlist(i)
               amin(i)    = dble(ztot(i)) / amin(i)
               ntot(i)    = atot(i) - ztot(i)
            end do

            98778   i = i - 1

            close(isotopedatabase_fh)

            if (i < isotot) then
               print *,"Table isotopedatabase has only ",i,&
                     " isotopes, but you want to read ",isotot," as specified in"//&
                     "parameter nnn in ppn_physics.input. STOP."
               stop
            end if

            if ( isotot >  nsp ) then
               print *, 'isotot, nsp = ', isotot, nsp
               stop "error networkI: parameter nsp too small"
            end if
         end if

#ifndef PPN
         call broadcast(ztot)         ; call broadcast(atot)              ; call broadcast(ntot)
         call broadcast(amin)         ; call broadcast(considerisoinlist)
         call broadcast(amin)         ; call broadcast(isotot)
         call broadcast_ch_arr(eltot)
#endif

   end subroutine read_isotopedatabase



   subroutine zis_create_hash()

         ! ^_^ create a hash table for the lookup of the ppn isotope names in the zis
         !     vector. See ispe function and the utils module for more details

         zis_name_to_idx = create_hash_table(zis)

   end subroutine zis_create_hash



   subroutine nuc_data_init()

#ifndef PPN
         call nuc_data_broadcasts()
#endif

         call zis_create_hash()

         call create_az_to_zis_hash()

   end subroutine nuc_data_init



   subroutine create_az_to_zis_hash()

         ! ^_^ creates a hash table for looking up a ppn isotope name from the
         !     (a,z) values. The answer will not be necessarily correct for the
         !     isomers, but since this will be used for masking reaclib rates
         !     that we are not interested in, it is not a problem because it will
         !     be inclusive of isomers. See also: function epsi

         integer, dimension(nsp) :: a_int, z_int
         integer :: i

         do i = 1, nsp
            a_int(i) = int(anetw(i))
            z_int(i) = int(znetw(i))
         end do

         az_to_zis_hash_idx = create_hash_table(a_int, z_int, zis)

   end subroutine create_az_to_zis_hash




   integer function ispe( name )
         ! ^_^ given a ppn isotope name, e.g. "U 235", returns the index of this isotope
         !     in the zis array.

         integer                 :: i, idx
         character(len=5)        :: name

         ispe = -1

         if (iachar(name(1:1)) > 91) then ! name is in lower case
            do i = 1, 2
               if (iachar(name(i:i)) >= 97 .and. iachar(name(i:i)) <= 122) then
                  name(i:i) = achar(iachar(name(i:i)) - 32)
               end if
            end do
         end if

         ! *** translations necessary to map isotope names from SE abundance
         !     initialisation to PP
         ! ^_^
         select case(name)
         case("DEUT")
            name = "H   2"
         case("DEUTR")
            name = "H   2"
         case("AL26G")
            name = "AL 26"
         case("GAMMA")
            name = "OOOOO"
         case("H   1")
            name = "PROT "
         case("PROTO")
            name = "PROT "
         case("NEUTR")
            name = "NEUT "
         case default
            continue
         end select

         ! ^_^ at the moment, the zis hash table is only written at the end of
         !     rnetw2007, but the ispe function is called within rnetw2007. I think that
         !     the has table should, if possible, be written at the beginning of the
         !     physics initialisation, i.e. start of rnetw2007, but for now I just
         !     revert to the old ispe method (non hash table, order n search) if the has
         !     table has not been written yet

         if (allocated(zis_name_to_idx% values)) then

            ! ^_^ hash lookup

            idx = hash(name)
            do i = 1, zis_name_to_idx% nvalsinbucket
               if (zis_name_to_idx% values(idx,i) == -1) then
                  ispe = -1
                  exit
               else if (name == zis(zis_name_to_idx% values(idx,i))) then
                  ispe = zis_name_to_idx% values(idx,i)
                  exit
               end if
            end do

         else

            ! ^_^ O(n) search

            do i = 1, nsp
               if (name == zis(i)) then
                  ispe = i
                  exit
               endif
            end do

         end if

   end function ispe



   function epsi(a, z)

         ! ^_^ inverse of ispe, i.e. it gives the ppn isotope name given A and Z

         integer :: a, z
         integer :: idx, i, ispe_index
         character(len=5) :: epsi
         character(len=2), parameter :: &
            D(0:112) = [ 'nn', 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', &
            'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
            'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', &
            'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', &
            'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', &
            'SN', 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', &
            'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', &
            'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', &
            'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', &
            'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', &
            'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', &
            'DS', 'RG', 'CN' ]

         if (allocated(az_to_zis_hash_idx% values)) then

            idx = Cantor_hash(a,z)

            do i = 1, az_to_zis_hash_idx% nvalsinbucket
               epsi = az_to_zis_hash_idx% values(idx, i)
               ispe_index = ispe(epsi)
               if (ispe_index == -1) then
                  !print *, 'species not found in epsi function', a, z
                  epsi = 'nojoy'
               else if (a == int(anetw(ispe_index)) .and. z == int(znetw(ispe_index))) then
                  exit
               else
                  epsi = 'nojoy'
               end if
            end do
         else 
            ! create name and return it
            if (a == 1 .and. z == 0) then
               epsi = "NEUT "
            else if (a == 1 .and. z == 1) then
               epsi = "PROT "
            else
               epsi = '     '
               write(epsi(:2),"(a2)") D(z)
               write(epsi(3:),"(i3)") a
            end if
         end if

   end function epsi



   !> knowing the element, this subr. provides the corrispondent atomic number
   subroutine atomic(nre1,netelt,netelf,netzt,netzf,ireac)
         implicit none
         integer ireac,nre1,i
         integer netzt(nre1),netzf(nre1)
         character(len=2) netelt(nre1),netelf(nre1)
         character(len=2), parameter :: &
               D(86) = [ 'nn', 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', &
               'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
               'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', &
               'ZN', 'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', &
               'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', &
               'SN', 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', &
               'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', &
               'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', &
               'HG', 'TL', 'PB', 'BI', 'PO', 'AT' ]

         do i=1,86
            if ( netelt(ireac)  ==  D(i) ) then
               netzt(ireac) = i-1
               exit
            end if
         end do

         do i=1,86
            if ( netelf(ireac)  ==  D(i) ) then
               netzf(ireac) = i-1
               exit
            end if
         end do
   end subroutine atomic


#ifndef PPN
   subroutine nuc_data_broadcasts()
         use communication

         ! ^_^ characters
         call broadcast_ch_arr(zis)

         !^_^ logicals
         call broadcast(considerisotope)
         call broadcast(considerreaction)

   end subroutine nuc_data_broadcasts
#endif


end module nuc_data
