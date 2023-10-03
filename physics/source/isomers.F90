!> @autor Rene Reifarth, Christoph Köppchen
!
!> This module calculate the effective decay constant between the groundstate and the isomere and
!> between the different isomeric states. It also computes the effective beta decay constant using
!> the given energy level files. These calculations determines the isomeric states and store the 
!> involving reactions and their deacy constants as a function of the temperature in the rates file. 
!> The isomers will be added to the network. All reactions even isomer to isomer reactions will be 
!> added to the network. During the runtime all the reaction rates for the isomers will be determined.
!
module isomers
   use array_sizes, only: nsp, nre
   use utils, only: r8, i4
   use constants
   use physics_knobs, only: isomer_choice
   implicit none

   private

   ! file input
   real(r8), parameter :: kev_to_t9 = 11604525e-9_r8

   integer :: isomer_filenumber, filenumber_rates, filenumber_selcect_isomers, &
         filenumber_decay_file, filenumber_element, filenumber_levelname
   character(len = 256) :: isomer_filename = '../NPDATA/isomerfiles/isomer+level'
   character(len = 256) :: filename_rates = '../NPDATA/isomerfiles/rates'
   character(len = 256) :: filename_select_isomers = '../NPDATA/isomerfiles/select_isomers'
   character(len = 256) :: decay_file
   character(len = 5) :: zis_a, zis_b, zis_gr, zis_gr1
   character(len = 3) :: arrow
   character(len = 1) :: isolevelchar
   character(len = 1), dimension(53) :: levelname

   real(r8) :: df

   ! variables for calculation of the effectiv decay rates and transitions

   character(len = 20)  element_names(0:200)

   real(r8) :: time_cut, lambda_cut, kt

   real(r8), dimension (:), allocatable, target :: &
         e_product, & ! Energy level in product nucleus
         e_level, &   ! Energy level of state (keV)
         g_level, &   ! degeneracy of level
         j_level, &   ! angular momentum J of level (h_bar)
         p_level, &   ! parity of level (+-1)
         lambda, &    ! terrestrial decay constant of state (1/s)
         temp_reac    ! temparature for all new beta decays and transitions

   real(r8), dimension (:,:), allocatable, target :: &
         transition, &                  ! terrestrial probability of internal transition as a function of initial level and final level
         decay, &                       ! terrestrial probability of decay as a function of initial level and reaction channel
         lambda_transition, &           ! terrestrial internal transition rate as a function of initial level and final level
         lambda_transition_wk, &        ! Weisskopf prediction of terrestrial internal transition rate as a function of initial level and final level
         lambda_decay, &                ! terrestrial decay rate as a function of initial level and reaction channel
         lambda_kt, &                   ! thermal excitation rate as a function of initial level and final level
         lambda_transition_effective, & ! effective stellar internal transition rate as a function of initial and final long-lived ( !  ! ) level
         lambda_decay_effective, &      ! effective stellar decay rate as a function of initial long-lived level and reaction channel
         all_reac                       ! store all reaction rates for beta decay and transitions
   integer(i4) :: &
         num_isomers, num_active_isomers, &
         number_reaction_types, number_levels, number_isomers, number_intermediate, &
         az, aa, level1, zz, za, level2, &
         number_new_reactions, number_of_all_reactions, &
         a, z, n, &                      ! A, Z, N of isotope under investigation
         internal_level_id(0:1000000), & ! gives internal level-id as a function of the external ID (external ID must be inside [0,1000000])
         internal_reaction_id(100)       ! gives internal reaction id as a function of the external ID (external ID must be inside [1,100])

   integer(i4), dimension(nre) :: &
         indi, &                       ! points to index of isomer-counterpart reaction from ground state reaction index
         flag_isomer, &                ! flag for isomer involvement:
                                       !     0: no isomer involvement
                                       !     1: isomer involved; this reaction involves the ground state
                                       !     2: isomer involved; this reaction involves the isomer
         inn, &                        ! points to the isomer (1-num_isomers) involved in this reaction
         new_reac                      ! flag for new reactions for beta decays and transitions

   integer, dimension (:), allocatable, target :: &
         z_product, &                   ! Charge of reaction product nucleus
         n_product                      ! Neutron number of reaction product nucleus



   public isomers_init, num_isomers, isomers_calculate_rates, num_active_isomers, isolevelchar, df
#ifndef PPN
   public isomers_broadcasts
#endif

contains

   !---
   !> @brief Before the runtime isomeric states and decay constants are being calculated. 
   !> Isomers and necessary reactions will be added. 
   !--- 
   subroutine isomers_init(an, zn, nums, numr, lab, labb, ilabb)
         use nuc_data
         use rates
         implicit none

         integer :: i, j, l, k, ki, jj, ll, &
               nums, &     ! number of species active in network
               numr, &     ! number of reactions acive in network without isomer reactions
               na          ! index of isomer in global species array
         integer :: ilabb(nre)
         real(r8) :: an(nsp), zn(nsp)
         character(len=5) :: spe, lab(nre), labb(nre)
         integer :: number_state, number_state1
         character(len=5), dimension(:), allocatable :: zis_mer, zis_mer1
         real(r8), dimension(:), allocatable :: isomer_spin, isomer_level_energy
         character(len=5) :: zis_gr, zis_gr1
         integer :: a_iso, z_iso, a_iso1, z_iso1


         if (isomer_choice == 0) then
            num_isomers = 0
            num_active_isomers = 0
            return
         end if

         call init_isomer_rates()

         open(newunit=isomer_filenumber,file = isomer_filename)
         read(isomer_filenumber,*)
         read(isomer_filenumber,*) num_isomers
         read(isomer_filenumber,*)
         close(isomer_filenumber)

         na = 0; ll = 0

         do i = 1, num_isomers !include all isomeric states
            call isomer_read(i, zis_gr, number_state, zis_mer, a_iso, z_iso)
            do j = 1, nums
               if (.not. considerisotope(j)) cycle
               do k = 1, number_state
                  if (zis_gr == zis(j)) then
                     spe      = zis_mer(k)
                     na       = nums + ll + k
                     zis(na)  = spe
                     an(na)   = dble(a_iso)   
                     zn(na)   = dble(z_iso)
                     considerisotope(na) = .true.
                  end if
               end do
            end do
            ll = ll + number_state
         end do

         ! update number of isomers (because network may not be large enough to include /all/ isomers)
         num_active_isomers = na - nums
         ! update number species to include isomers
         nums = na

         flag_isomer(:) = 0; inn(:) = 0; indi(:) = 0; jj = 0; na = 0; ll = 0

         do j = 1, num_isomers !include reaction rates from nuklid to isomer and vice versa
            call isomer_read(j, zis_gr, number_state, zis_mer, a_iso, z_iso)
            do l = 1, numr
               if (.not. considerreaction(l)) cycle
               do k = 1, number_state
                  if (zis(k1(l)) /= zis_gr .and. zis(k7(l)) /= zis_gr) cycle
                  jj          = jj + 1
                  spe         = zis_mer(k)
                  na          = numr + jj
                  k2(na)      = k2(l)
                  k3(na)      = k3(l)
                  k4(na)      = k4(l)
                  k5(na)      = k5(l)
                  k6(na)      = k6(l)
                  k8(na)      = k8(l)
                  considerreaction(na) = considerreaction(l)
                  ilabb(na)   = ilabb(l)
                  lab(na)     = lab(l)
                  labb(na)    = labb(l)
                  indi(l)     = na
                  flag_isomer(l)  = 1
                  flag_isomer(na) = 2 
                  inn(na)      = l  
                  if (zis(k1(l)) == zis_gr) then
                     k1(na) = ispe(spe)
                     k7(na) = k7(l)
                     if (k1(l) == k7(l)) k7(na) = k1(na)
                  else if (k1(l) /= k7(l)) then
                     k1(na) = k1(l)
                     k7(na) = ispe(spe)
                  end if
               end do
            end do
         end do


         do j = 1, num_isomers !include reaction rates from isomer to isomer
            call isomer_read(j, zis_gr, number_state, zis_mer, a_iso, z_iso)
            do i = 1, num_isomers
               call isomer_read(i, zis_gr1, number_state1, zis_mer1, a_iso1, z_iso1)
               if (i /= j) then
                  do l = numr, numr+jj
                     if (.not. considerreaction(l)) cycle
                     do k = 1, number_state
                        if (zis(k1(l)) == zis_mer(k) .and. zis(k7(l)) == zis_gr1)then
                           do ki = 1, number_state1
                              jj          = jj + 1
                              spe         = zis_mer1(ki)
                              na          = numr + jj
                              k2(na)      = k2(l)
                              k3(na)      = k3(l)
                              k4(na)      = k4(l)
                              k5(na)      = k5(l)
                              k6(na)      = k6(l)
                              k8(na)      = k8(l)
                              considerreaction(na) = considerreaction(l)
                              ilabb(na)   = ilabb(l)
                              lab(na)     = lab(l)
                              labb(na)    = labb(l)
                              flag_isomer(na)      = 2 
                              inn(na)      = l
                              k1(na) = k1(l)
                              k7(na) = ispe(spe)
                           end do
                        end if
                     end do
                  end do
               end if
            end do
         end do

         call identify_number_reaction(number_new_reactions)

         open(newunit = filenumber_rates, file=filename_rates) 

         read(filenumber_rates,*)

         do l = 1, number_new_reactions ! include transitions
            read(filenumber_rates,'(2(1x,I3),1x,I2,2(1x,I3),1x,I2,1x,A5,1x,A3,1X,A5)')&
                  & az, aa, level1, zz, za, level2, zis_a, arrow, zis_b

            if (level2 /= 0) then   
               na = na + 1
               k1(na) = ispe(zis_a)
               k2(na) = 1
               k3(na) = ispe("OOOOO")
               k4(na) = 0
               k5(na) = ispe("OOOOO")
               k6(na) = 0
               k7(na) = ispe(zis_b)
               k8(na) = 1
               considerreaction(na) = .true.
               ilabb(na)   = 22
               lab(na)     = "?????"
               labb(na)    = "(+,+)"
               flag_isomer(na) = 2
               new_reac(na) = l
            else
               do i = 1, na
                  if(zis(k1(i))==zis_a .and. zis(k7(i))==zis_b .and. zis(k3(i))=='OOOOO' .and. zis(k5(i))=='OOOOO')then
                     new_reac(i) = l
                  end if
               end do    
            end if

            do i = 1, 1000
               read(filenumber_rates,*) 
            end do
         end do

         close(filenumber_rates)
    
         allocate(temp_reac(1000), all_reac(number_new_reactions, 1000))

         open(newunit = filenumber_rates, file=filename_rates)
         read(filenumber_rates,*)
         do l = 1, number_new_reactions
            read(filenumber_rates,*)
               do i = 1, 1000
                  read(filenumber_rates, '(2(es23.15))') temp_reac(i), all_reac(l,i)
               end do   
         end do
         close(filenumber_rates)

         temp_reac = temp_reac*kev_to_t9

         number_of_all_reactions = na
         numr = na

   end subroutine isomers_init

  
   !---
   !> @brief This subroutine picks out the name of a choosen nuclide as well as its isomere, the count 
   !> of the isomers, mass number and charge number of the file isomer+level
   !---
   subroutine isomer_read(chosen_nuclid, zis_gr, number_state, zis_mer, a_iso, z_iso)
         implicit none
         integer :: chosen_nuclid, number_state, number_isomer
         integer :: a_iso, z_iso
         character(len = 5) :: zis_gr, dummy
         character(len = 5), dimension(:), allocatable :: zis_mer
         real(r8), dimension(:), allocatable :: isomer_spin, isomer_level_energy
         integer :: i, ii, k

         open(newunit = isomer_filenumber,file=isomer_filename)
         read(isomer_filenumber,*)
         read(isomer_filenumber,*) number_isomer
         read(isomer_filenumber,*)

         do ii = 1, chosen_nuclid-1
            read(isomer_filenumber,'(1x,A5,2x,I2)') dummy, number_state
            do k = 1, number_state
               read(isomer_filenumber,*)   
            end do
         end do

         read(isomer_filenumber,'(1x,A5,2x,I2,2x,I3,2x,I3)') zis_gr, number_state, z_iso, a_iso 

         if (allocated(zis_mer)) deallocate(zis_mer)

         allocate(zis_mer(number_state))

         do i = 1, number_state
            read(isomer_filenumber,'(1x,A5,1x,es23.15,1x,F7.2)') zis_mer(i)
         end do

         close(isomer_filenumber)

   end subroutine isomer_read


   !---
   !> @brief The subrotine identify the number of transitions and new beta decays of the file rates
   !---
   subroutine identify_number_reaction(number_reaction)
         implicit none
         integer :: i, number_reaction, iostats = 0

         number_reaction = 0   

         open(newunit = filenumber_rates, file = filename_rates)

         do while (iostats == 0)
            read(filenumber_rates,*, iostat = iostats)
            do i = 1, 1000
               read(filenumber_rates,*,iostat=iostats)
            end do
            number_reaction = number_reaction + 1
         end do 

         number_reaction = number_reaction - 1

         close(filenumber_rates)

   end subroutine identify_number_reaction

   !---
   !> @brief Interpolation between two decay constants for a given temperature.
   !---
   subroutine interpolate_rate(t9, v_new, a, b)
         implicit none
         real(r8) :: v_new, t9
         integer :: i, t_location
         real(r8) :: m, b1
         real(r8), dimension(1000) :: a, b

         t_location = 0

         do i = 1, 1000
            if(a(i) <= t9)then
               t_location = t_location +1
            end if
         end do 

         m = (b(t_location+1) - b(t_location))/((a(t_location + 1) - a(t_location)))

         b1 = b(t_location) - m*a(t_location)

         v_new = m * t9 + b1

   end subroutine interpolate_rate


   !---
   !> @brief During the runtime reaction rates are set for a given temperature.
   !--- 
   subroutine isomers_calculate_rates(t9, numr, lab, labb, ilabb)
         use rates, only: v, k1, k7, k3, k5
         use nuc_data
         implicit none

         integer :: j, l, numr, r_typ, idx, i
         integer :: na ! index of isomer reaction in master reaction list
         real(r8) :: rgdummy, ridummy, t9, v_new
         integer :: ilabb(nre)
         character(len=5) :: lab(nre), labb(nre)


         do l = 1, number_of_all_reactions
            if(flag_isomer(l)==2)then  
               na = inn(l)
               v(l) = v(na)
            end if
         end do

         do l = 1, number_of_all_reactions
            do j = 1, number_new_reactions
               if(new_reac(l) == j)then
                  call interpolate_rate(t9, v_new, temp_reac, all_reac(j,:))
                  v(l) = v_new
               end if
            end do
         end do

   end subroutine isomers_calculate_rates


   !---
   !> @brief This subroutine creates the files isomer+level and rates for the nuclides of the file 
   !> select_isomers. The file isomer+level consists of the name of the nuclide, name of the isomers,
   !> number of the isomers, mass number and charge number. The file rates consists of reaction, decay 
   !> constants as a function of the temperature.
   !---
   subroutine init_isomer_rates()
         implicit none
         integer :: i, number_of_isomer, kt_int, j, k, zero, number_of_levels
         character(len = 5), dimension(:), allocatable :: zis
         integer, dimension(:), allocatable :: az, aa
         character(len = 5) :: testchar, dummy1, dummy2
         character(len = 3) :: eledummy

         open(newunit = filenumber_levelname, file = '../NPDATA/isomerfiles/levelname')
    
         do i = 1, 53
         read(filenumber_levelname,'(A1)') levelname(i)
         end do

         close(filenumber_levelname)  

         open(newunit = filenumber_selcect_isomers, file=filename_select_isomers, status="old")

         do i = 1, 3
            read(filenumber_selcect_isomers, *) !Header
         end do

         read(filenumber_selcect_isomers,*) number_of_isomer

         allocate(az(number_of_isomer),aa(number_of_isomer),zis(number_of_isomer))

         do i = 1, number_of_isomer
            read(filenumber_selcect_isomers,'(1x,A5,2x,I3,2x,I3)') zis(i), aa(i), az(i)
         end do

         close(filenumber_selcect_isomers)

         open(newunit = filenumber_rates,file=filename_rates)

         write(filenumber_rates,*) "   az   aa   level   zz   za   level   reaction"  

         open(newunit = isomer_filenumber,file=isomer_filename)
         write(isomer_filenumber,*) "number of nuklids considered as isomer"
         write(isomer_filenumber,*) number_of_isomer
         write(isomer_filenumber,*) "isotope    number of iosomer + ground state"

         do k = 1,number_of_isomer !Loop over all isotope

            call initfiles(zis(k))
            call initdecay 

            number_of_levels = 0

            do i=1, number_levels
               if (lambda(i) < lambda_cut) then  
                  number_of_levels = number_of_levels + 1
               end if
            end do 

            write(isomer_filenumber,'(1x,A5,2x,I2,2x,I3,2x,I3)') zis(k), number_of_levels - 1, az(k), aa(k)  
            ! call calculate_rates
            do i=1, number_levels
               if (lambda(i) < lambda_cut) then            ! isomer/gs
                  do j=1, number_reaction_types 
                     testchar = zis(k)
                     dummy1 = testchar(1:2) // isolevelchar(i) // testchar(4:5)
                     eledummy = "   "
                     write(eledummy,"(I3)") (n_product(j) + z_product(j))
                     dummy2 = trim(element_names(z_product(j))) // eledummy
                     zero = 0
                     write(filenumber_rates,'(2(1x,I3),1x,I2,2(1x,I3),1x,I2,1x,A5,1x,A3,1X,A5)') &
                           & z, a, i, z_product(j), (z_product(j)+n_product(j)), zero , dummy1, "-->", dummy2       
                     !write(131,*) "Decay from Level", i, "   ", zis(k)
                     do kt_int=1,1000
                        kt = kt_int
                        call calculate_rates
                        write(filenumber_rates, '(2(es23.15))') kT, lambda_decay_effective(i,j)
                     end do 
                  end do

                  do j=1, number_levels     
                     if (lambda(j) < lambda_cut .and. i /= j) then            ! isomer/gs
                        testchar = zis(k)
                        dummy1 = testchar(1:2) // isolevelchar(i) // testchar(4:5)
                        dummy2 = testchar(1:2) // isolevelchar(j) // testchar(4:5)
                        write(filenumber_rates,'(2(1x,I3),1x,I2,2(1x,I3),1x,I2,1x,A5,1x,A3,1X,A5)') &
                              &z, a, i, z, a, j, dummy1, "-->",dummy2 
                        !write(131,*) "transition from ", i, "to", j 
                        do kt_int=1,1000
                           kt = kt_int
                           call calculate_rates
                           write(filenumber_rates, '(2(es23.15))') kT, lambda_transition_effective(i,j)
                        end do 
                     end if  
                  end do

                  if(i /= 1)then
                     write(isomer_filenumber,'(1x,A5,1x,es23.15,1x,F7.2)') dummy1
                  end if
               end if  
            end do

            deallocate (z_product, n_product, e_product, g_level, e_level, j_level, p_level)
            deallocate (lambda, transition, decay, lambda_transition_wk, lambda_transition, lambda_transition_effective)
            deallocate (lambda_kt, lambda_decay, lambda_decay_effective)

         end do

         close(filenumber_rates)
         close(isomer_filenumber)


   end subroutine init_isomer_rates


   !---
   !@brief This subroutine reads the file element_symbols.txt and stores in element_names.
   !---
   subroutine initfiles(input)

         implicit none

         character(len=5) :: input
         integer :: iostats, i, j

         decay_file = "../NPDATA/isomerfiles/energylevel_files/" // input//".sch"

         open(newunit = filenumber_decay_file,file=decay_file,status='unknown')
         read(filenumber_decay_file,*) i,j,time_cut
         !write(*,*) 'Isomer time cut    : ', time_cut, ' s'
         close(filenumber_decay_file)

         lambda_cut = dlog(2.d0)/time_cut
         open(newunit = filenumber_element,file='../NPDATA/isomerfiles/element_symbols.txt',status='old')
         iostats=0
         i = 0
         element_names(0) = 'NN'

         do while (iostats == 0)
            i=i+1
            read(filenumber_element,*,iostat=iostats) element_names(i)
         end do 

         close(filenumber_element)

   end subroutine initfiles



   subroutine initdecay
         implicit none

         integer :: level

         internal_level_id    = -1
         internal_reaction_id = -1
         call initlevels
         allocate (lambda_transition(number_levels,number_levels))
         allocate (lambda_transition_effective(number_levels,number_levels))
         allocate (lambda_kt(number_levels,number_levels))
         allocate (lambda_decay(number_levels,number_reaction_types))
         allocate (lambda_decay_effective(number_levels,number_reaction_types))

         lambda_transition = 0._r8
         lambda_transition_effective = 0._r8
         lambda_kt = 0._r8
         lambda_decay = 0._r8
         lambda_decay_effective = 0._r8

         call init_wk

         do level=1,number_levels
            lambda_transition(level,:) = lambda(level) * transition(level,:)
            lambda_decay(level,:)      = lambda(level) * decay(level,:)
         end do

         call get_rates_wk

         do level=1,number_levels
            if (lambda(level) <= lambda_cut) then
               number_isomers      = number_isomers + 1
               print*,'Level ',level,' - isomer, half-life time: ',0.693/lambda(level),' s'
            else
               number_intermediate = number_intermediate + 1
               !print*,'Level ',level,' - intermediate state, half-life time: ',0.693/lambda(level),' s'
            end if  
         end do

         !print*,'Number of long-lived states  :', number_isomers
         !print*,'Number of intermediate states:', number_intermediate
         ! stop

   end subroutine initdecay


   !---
   !> @brief This subroutine reads the energy level file for the calculations of the effective deacay constants
   !---
   subroutine initlevels

         implicit none
         integer :: i, j, reaction_id, dz, dn, level_id
         real(r8) :: level_energy, level_j, t12, branching, norm, level_p
         integer :: nbr_target_levels

         open(newunit = filenumber_decay_file,file=decay_file,status='old')
         read(filenumber_decay_file,*) z, n
         a = z+n
         !print *,'A, Z, N: ', a, z, n 
         read(filenumber_decay_file,*) number_reaction_types
         !print *,'Number of reaction types: ', number_reaction_types

         if (number_reaction_types > 0) then
            allocate (z_product(number_reaction_types))
            allocate (n_product(number_reaction_types))
            allocate (e_product(number_reaction_types))

            do i=1,number_reaction_types
               read(filenumber_decay_file,*) reaction_id, dz, dn, level_energy
               if (reaction_id == 0) then
                  print*,'Reaction-ID must not be 0 - this is reserved for internal transition. Check input.'
                  stop
               end if  
               if (internal_reaction_id(reaction_id) == -1) then
                  internal_reaction_id(reaction_id) = i
               else
                  print*,'Reaction-ID ',reaction_id, ' not unique. Check input.'
                  stop
               end if  
               z_product(i) = z + dz
               n_product(i) = n + dn
               e_product(i) = level_energy
               !print*,reaction_id, z_product(internal_reaction_id(reaction_id)), n_product(internal_reaction_id(reaction_id)), &
               !level_energy,internal_reaction_id(reaction_id)
            end do
            read(filenumber_decay_file,*) number_levels
            !print*,'Number of levels:', number_levels
            allocate (e_level(number_levels))
            allocate (g_level(number_levels))
            allocate (j_level(number_levels))
            allocate (p_level(number_levels))
            allocate (lambda(number_levels))
            allocate (transition(number_levels,number_levels))
            allocate (decay(number_levels,number_reaction_types))
            allocate (lambda_transition_wk(number_levels,number_levels))

            lambda_transition_wk = 0._r8
            transition = 0._r8
            decay = 0._r8

            do i=1,number_levels
               read(filenumber_decay_file,*) level_id, level_energy, level_j, level_p, t12, nbr_target_levels
               !write(*,*) level_id
               if (internal_level_id(level_id) == -1) then
                  internal_level_id(level_id) = i
               else
                  print*,'Level-ID ',level_id, ' not unique. Check input.'
                  stop
               end if  
               e_level(i) = level_energy
               j_level(i) = level_j
               p_level(i) = level_p
               g_level(i) = 2._r8*j_level(level_id) + 1._r8
               if (t12 <= 0._r8) then
                  lambda(i)  = -1.d0
               else  
                  lambda(i)  = dlog(2._r8) / t12
               end if  
               norm = 0._r8
               do j=1, nbr_target_levels
                  read(filenumber_decay_file,*) reaction_id, level_id, branching
                  if (nbr_target_levels == 1) branching=1._r8
                  if (reaction_id == 0) then  !  IT
                     if (internal_level_id(level_id) > -1) transition(i,internal_level_id(level_id)) = branching
                  else if  (reaction_id == 1) then  !  Reaction changing N,Z
                     if (internal_reaction_id(level_id) == -1) then
                        print*,'Reaction-ID ',level_id , ' not defined. Check input.'
                        stop
                     end if  
                     decay(i,internal_reaction_id(level_id)) = branching
                  else
                     print*,'something wrong with input. Type can only be 0 or 1, not ',  reaction_id
                     stop
                  end if
                  norm = norm + branching
               end do
               if (norm == 0._r8 .and. nbr_target_levels > 0) then
                  print*,'Something wrong with branchings - sum equals zero for level:', i
               else 
                  if (norm > 0._r8) then
                     transition(i,:) = transition(i,:) / norm
                     decay(i,:)      = decay(i,:) / norm
                  end if  
               end if  
            end do
         end if

         close(filenumber_decay_file)
         do i=2, number_levels
            if (e_level(i) <= e_level(i-1)) then
               print*,'Please sort levels according to energy.'
               stop
            end if  
         end do

   end subroutine initlevels


   !---
   !> @brief Initialisation of Weisskopf approximation
   !---
   subroutine init_wk
         implicit none

         integer :: i, j
         real(r8) :: norm

         !print*,'was here --------------------init_wk'
         do i=1,number_levels
            do j=i-1,1,-1
               call get_weisskopf(i,j)
            end do
            if (lambda(i) == -1._r8) then
               !print*,'Will assume Weisskopf approximation for half-live of state', i
               lambda(i) = sum(lambda_transition_wk(i,:))
               !print*,'WK - half life time:', dlog(2.d0)/lambda(i), ' s'
            end if
         end do

   end subroutine init_wk


   !---
   !> @brief This subroutine gives the rate for a transition using Weisskopf approximation.
   !---
   subroutine get_rates_wk
         implicit none

         integer :: i, j
         real(r8) :: norm

         norm = 0._r8
         do i=1,number_levels
            do j=i-1,1,-1
               if (transition(i,j) == 0._r8) then 
                  lambda_transition(i,j) = lambda_transition_wk(i,j)
                  !print*,'Will assume Weisskopf approximation for transition of state', i, ' to state',j
                  !print*,'WK - half life time:', dlog(2.d0)/lambda_transition_wk(i,j), ' s'
               end if  
            end do
            norm = sum(lambda_transition(i,:)) + sum(lambda_decay(i,:))
            transition(i,:) = lambda_transition(i,:)/norm
         end do

   end subroutine get_rates_wk


   !---
   !> @brief Calculation of the transition of two states i and j by Weisskopf approximation
   !---
   subroutine get_weisskopf(i,j)
         implicit none

         integer :: i,j, reaction_id, dz, dn, level_id, ll, hl
         real(r8) :: energy,l, h, c, alpha, mp, wk, hwz, r, r0
         logical :: l_even, dp, ttype

         ll = min(i,j)
         hl = max(i,j)
         h     = 6.58211e-22_r8                     ! MeV s
         c     = 299792458._r8                    ! m/s
         alpha = 1._r8/137._r8                         ! cgs
         mp    = 938.3_r8/c**2                      ! MeV*s^2/m^2
         r0    = 1.26_r8                         ! fm
         r = r0*(a**(1._r8/3._r8))/1e15_r8           ! m
         energy = ( e_level(hl) - e_level(ll) )  / 1.e3_r8 
         l      = dabs(j_level(hl) - j_level(ll))
         if (l > 0) then
            !print*,'Transition energy:', energy, ' MeV'
            l_even = (modulo(l,2._r8) == 0._r8)
            dp     = (p_level(ll)*p_level(hl) > 0._r8)                  ! dp = true, if no parity change, dp = false if parity change
            ttype  = (l_even .and. dp ) .or. (.not.l_even .and. .not.dp)   ! dp = true if E-transition, dp = false if M-transition
            if (ttype) then
               !print*,'Transition multipolarity:  E', int(l)
               wk = h*((2._r8*(l+1._r8))/(l*(df(2._r8*l+1._r8))**2))*(9._r8/(l+3._r8)**2)*((energy/(h*c))**(2._r8*l+1._r8))* &
                     alpha*c*(R**(2._r8*l)) ! MeV
            else 
               !print*,'Transition multipolarity:  M', int(l)
               wk = h*((20._r8*(l+1._r8))/(l*(df(2._r8*l+1._r8))**2))*(9._r8/(l+3._r8)**2)*((energy/(h*c))**(2._r8*l+1._r8))* &
                     alpha*c*(R**(2._r8*l-2._r8))*((h/(mp*c))**2) ! MeV
            end if
            !print*,'Weisskopf unit: ', wk,' MeV'
            !hwz = dlog(2.d0)*h/wk
            !print*,'Half life time: ', hwz, ' s'
            lambda_transition_wk(i,j) = wk/h
         else
            lambda_transition_wk(i,j) = 0._r8
         end if

   end subroutine get_weisskopf



   subroutine calculate_rates
         implicit none

         integer :: i, j, k, l
         real(r8) :: ll, sum_k, sum_l

         call calculate_thermal_excitation 
         ! effective transitions between long-lived states
         lambda_transition_effective = 0._r8
         do i=1, number_levels
            if (lambda(i) < lambda_cut) then            ! isomer/gs
               do j=1, number_levels                      
                  if (lambda(j) < lambda_cut .and. j /= i) then        ! isomer/gs
                     ll = lambda_kT(i,j) + lambda_transition(i,j)
                     sum_l=0._r8
                     do l=1, number_levels
                        if (lambda(l) >= lambda_cut) then            ! intermediate state
                           sum_k=0._r8
                           do k=1,l-1
                              if (lambda(k) >= lambda_cut) then            ! intermediate state
                                 sum_k = sum_k + transition(l,k)*transition(k,j)
                              end if
                           end do
                           sum_l = sum_l + lambda_kt(i,l)*(transition(l,j)+sum_k) 
                        end if 
                     end do
                     lambda_transition_effective(i,j) = ll+sum_l
                  end if
               end do
            end if
         end do  

         !effective decays of long-lived states

         do i=1, number_levels
            if (lambda(i) < lambda_cut) then            ! isomer/gs
               do j=1, number_reaction_types    
                  ll = lambda_decay(i,j)
                  sum_l=0._r8
                  do l=1, number_levels
                     if (lambda(l) >= lambda_cut) then            ! intermediate state
                        sum_k=0._r8
                        do k=1,l-1
                           if (lambda(k) >= lambda_cut) then            ! intermediate state
                              sum_k = sum_k + transition(l,k)*decay(k,j)
                           end if
                        end do
                        sum_l = sum_l + lambda_kt(i,l)*(decay(l,j)+sum_k)

                     end if 
                  end do
                  lambda_decay_effective(i,j) = ll+sum_l

               end do
            end if
         end do 


   end subroutine calculate_rates



   subroutine calculate_thermal_excitation
         implicit none

         integer :: i, j

         lambda_kt = 0._r8
         do i=1,number_levels
            do j=1,number_levels
               if (lambda_transition(i,j) > 0._r8) then
                  lambda_kt(j,i) = lambda_transition(i,j) * g_level(i)/g_level(j)
                  lambda_kt(j,i) = lambda_kt(j,i) * dexp((e_level(j)-e_level(i))/kt)
               end if  
            end do
         end do
   end subroutine calculate_thermal_excitation



#ifndef PPN
   subroutine isomers_broadcasts()
         use communication
         implicit none

         call broadcast(indi)
         call broadcast(flag_isomer)
         call broadcast(inn)

   end subroutine isomers_broadcasts
#endif


   !---
   !> @brief Support function for the naming of the isomers
   !---
   function isolevelchar(level)
         implicit none

         integer :: level

         if(level == 1 .and. a < 100)then
            isolevelchar = levelname(1)
         elseif(level == 1 .and. a >=100 .and. a<200)then
            isolevelchar = "1"
         elseif(level == 1 .and. a >200)then
            isolevelchar = "2"
         else
            isolevelchar = levelname(level)
         end if

   end function isolevelchar



   recursive real(r8) function df(n) result(res)
         use utils, only: r8
         implicit none
         real(r8) n

         if (n < 2) then
            res = 1._r8
         else
            res = n*df(n-2._r8)
         end if

   end function df    


end module isomers
