!> BASEL reaclib comments:                                                                  
!>------------------------                                                                  
!> *** reaclib, Basel version (Rauscher & Thielemann, version available in www.kadonis.org, 
!>     15th April 2008)                                                                     
!> *** warning = there are not all the partition functions reqired for reaclib              
!> *** network. Is there only a mistake? Are available those info ?                         
!> *** 5412 + 1(mg18) in winvn.nosmo.correct (data of mg19)                                 
!> *** 5412 + 1(mg18) + 1(po277) <-- wrong! (data of po276)                                 
!> *** 5412 + 2(mg18,po277) + 1(at280) <-- wrong! (data of at279)                           
!> *** 5412 + 3(mg18,po277,At280) + 1(f15) <-- wrong! (data of f14)                         
!> *** 5416 + 11(f27-f37) <-- wrong! (data of f26)                                          
!>                                                                                          
!>
!> JINA reaclib comments:                                                                   
!>-----------------------                                                                   
!> *** warning = there are not all the partition functions reqired for reaclib              
!> *** network. Is there only a mistake? Are available those info ?                         
!> *** 5412 + 1(mg18) in winvn.nosmo.correct (data of mg19)                                 
!> *** 5412 + 1(mg18) + 1(po277) <-- wrong! (data of po276)                                 
!> *** 5412 + 2(mg18,po277) + 1(at280) <-- wrong! (data of at279)                           
!> *** 5412 + 3(mg18,po277,At280) + 1(f15) <-- wrong! (data of f14)                         
!> *** 5416 + 11(f27-f37) <-- wrong! (data of f26)                                          
!>
!> @todo during init, when masks are created, make new data arrays only contining the reaclib data
!> that we are actually going to need for the biggest network that we have, and deallocate
!> (if possible) the main data arrays. This will save memory consumption, encourage more contiguous memory
!> access and improve computation time. This is not a trivial re-write but will improve code speed
!> and hopefully, transparency
!> @todo set array sizes (either iAtdim... or iCfdim...) in init based in pIDX_RCLB only once,
!>       and then use this everywhere else, instead of having all the preprocessor statements I have
!>       put in now for #if pIDX_RCLB...

module reaclib
   use communication

   use array_sizes, only: ndim_rl, nnpartdim, partdim, iAtdim, i282dim, iCfdim, i325dim, nre
   use reaction_info, only: num_rtypes
   use utils, only: r8, wallclocktime, table_input, hash, fmtprofile, int_hash_table, &
         create_hash_table, irate, idrdt, idrdd
   use nuc_data, only: epsi, considerreaction
   use physics_knobs, only: index_reaclib, use_cache
   use constants

   implicit none
   private

   ! parameters
   integer, parameter   :: nchap = 11 ! number of chapters
   real(r8), parameter  :: t9_lowlimit_reaclib = 0.01_r8

   ! names of species for various things
   character(len=5), dimension(ndim_rl) :: e1, e2, e3, e4, e5, e6
   character(len=5) :: spe_part(nnpartdim+1), reaclib_name
   ! for hashing
   type(int_hash_table)  :: reaclib_name_to_idx

   ! reaclib reaction data: Q-value, the 7 coefficients, reference, rflags,
   !                        chapter, reactants' a, z, state
   real(r8), dimension(ndim_rl) :: qmev, a0, a1, a2, a3, a4, a5, a6
   character(len=4), dimension(ndim_rl)  :: rlreac
   character(len=1), dimension(ndim_rl), public  :: rflag2
   character(len=1), dimension(ndim_rl)  :: rflag1, chapter
   integer, dimension(ndim_rl) :: an1, an2, an3, an4, an5, an6, &
         zn1, zn2, zn3, zn4, zn5, zn6, &
         ista1, ista2, ista3, ista4, ista5, ista6
   logical, public :: is_reverse(ndim_rl)
   ! ^_^ for merging back into ppn network (?)
   !     nireac(ndim_rl): reaction type (1-15; see definitions in ppn_physics.F)
   integer, dimension(ndim_rl) :: ista, aan, zzn, istaf, aanf, zznf, nireac
   real(r8), dimension(ndim_rl), public :: rv
   integer                        :: jtotdim

   ! ^_^ give reaclib rate index (e.g. the index of a rate in r0 or reacv rate
   !     arrays, and get back the index of this reaction in the rv(:) rate vector
   !     that is used to merge reaclib rates into the ppn network's v(:) rate
   !     vector
   integer, dimension(ndim_rl)   :: map_reaclib_to_ppn

   ! ^_^ give ppn rate index (e.g. the index of a rate in rv(:) rate vector
   !     that is used to merge reaclib rates into the ppn network's v(:) rate
   !     vector, and get back the index of this reaction in the reaclib arrays,
   !     e.g. r0(:) or reacv(:)
   integer, dimension(ndim_rl), public   :: map_ppn_to_reaclib

   ! ^_^ and finally, the map from vbranch (dimension: nre) to complete the mapping
   !     back into the main ppn rates vector; i.e. the inverse of locate_reaclib_in_vbranch
   integer, dimension(ndim_rl), public :: locate_vbranch_in_reaclib_rv

   !     r0(:) contains the exponent a0(i) + a1(i)/t9 + a2(:)/t9^[...]
   !     reacv(:) is sum of all contributions to the rate exp(r0(a)) + exp(r0(b)) 
   !     pfactor is the partition function factor that the rate reacv is multiplied by
   real(r8), dimension(ndim_rl)  :: r0, reacv, pfactor, dr0dT, dreacvdT
   ! ^_^ reaclib_mask is a logical array that is true for reaclib reaction
   !     contributions(i.e. elements in the a0(:), a1(:), ..., a6(:), r0(:)
   !     arrays) that are to be used in the ppn network, based on the
   !     reactant a and z. reaclib_reacv_mask is the same but for complete
   !     reactions
   logical, dimension(ndim_rl) :: reaclib_r0_mask, reaclib_reacv_mask
   ! ^_^ pointer to tell us which rate (element in reacv array) a
   !     resonance (element in r0 array) contributes to
   integer, dimension(ndim_rl) :: resonance_to_rate

   ! taking into account that reaclib calculates in different steps resonances
   ! contributions and non-resonant rate, here below are the parent species,
   ! daughter species, reference (do I need rflags as well? I do not think for
   ! now) and Q-value rescaled for the right index (max number = ndim_rlt <
   ! ndim_rl). Those info below must be "saved" with costants for rate formula.
   integer :: ndim_rlt
   character(len=5), dimension(ndim_rl), public :: re1, re2, re3, re4, re5, re6
   character(len=4), dimension(ndim_rl) :: rrlreac
   character(len=1), dimension(ndim_rl) :: rchapter
   logical :: rrev(ndim_rl) !< if .true., reaction has flag "v", meaning it is a reverse rate

   ! reaclib nuclear data: A, Z, N, state, mass excess
   ! spinpart(i)          : ground state spin
   ! rt9(jj)              : temperature coordinates of partition function table (0.01 < T9 < 10)
   ! partnum(i,nnpartdim) : partition functions given in the grid temperature rt9
   ! partinter(nnpartdim) : is the interpolated partition functions for a given temperature
   !> WHY IS THIS nnpartdim+1, instead of nnpartdim???
   integer, dimension(nnpartdim+1)        :: anumpart, znumpart, istate
   integer, dimension(nnpartdim)          :: nnumpart
   real(r8), dimension(nnpartdim)         :: spinpart
   real(r8), dimension(nnpartdim,partdim) :: partnum
   real(r8), dimension(partdim)           :: rt9
   real(r8), dimension(nnpartdim)         :: exmass
   real(r8), dimension(nnpartdim+1)       :: partinter
   logical, dimension(nnpartdim+1)        :: reaclib_partfunc_mask
   ! ipart(i282dim+1,iAtdim+2,1) is for "empty specie" (+2 because neutron
   ! cannot have subscript Z = 0 in ipart, so Z + 1 as default subscript of
   ! ipart for all the species).
   
   ! this 4D array ntarget resembles reaclib_ntrans, which is used to map the reaclib
   ! rates back into the ppn network, but ntarget is calculated but never used
#if pIDX_RCLB == 3
   integer :: ipart(i325dim+1,0:iCfdim+1,3)
   integer :: reaclib_ntrans(i325dim,0:iCfdim,3,num_rtypes)
   integer :: ntarget(i325dim,0:iCfdim,3,num_rtypes)
#else
   integer :: ipart(i282dim+1,0:iAtdim+1,3)
   integer :: reaclib_ntrans(i282dim,0:iAtdim,3,num_rtypes)
   integer :: ntarget(i282dim,0:iAtdim,3,num_rtypes)
#endif

   ! dummy character variables for reading reaclib files
   character(len=1) :: cdummy
   character(len=5) :: cdummy1, cdummy2, cdummy3, cdummy4, cdummy5

   ! this seems to be a flag array that is 1 where we are interested in
   ! using the rate
   integer     :: mdummy(ndim_rl)

   ! caching/preprocessing
   logical     :: preprocessing = .false., cache_exists = .false.
   real(r8), allocatable :: reaclib_cache(:,:), drates(:)
   integer, parameter :: num_t9_preproc = 5001
   real(r8) :: dt9_preproc, t9_preproc(num_t9_preproc)


   ! timing
   real(r8)    :: tti1, tti2

   public reaclib_init, reaclib_wrapper, reaclib_interpolate_partition_functions, &
         reaclib_calculate_rates, index_reaclib, aan, zzn, &
         aanf, zznf, ista, istaf, reaclib_ntrans, nireac, jtotdim, &
         ndim_rlt, exmass, anumpart, znumpart, istate, spinpart, &
         partinter, partnum, reaclib_create_masks, reaclib_name, &
         reaclib_preprocessor, reacv

contains


   subroutine reaclib_create_hash()
         integer :: idx, i, j

         reaclib_name_to_idx = create_hash_table(spe_part)
   end subroutine reaclib_create_hash


   subroutine reaclib_wrapper ( t9in, rhoin, yein )
         use utils
         use physics_knobs, only: use_cache
         use array_sizes, only: ndim_rl
         implicit none
         real(r8), intent(in) :: t9in, rhoin, yein
         real(r8) :: t9, rho, ye
         ! pass
         t9 = t9in; rho = rhoin; ye = yein
         if ( .not. use_cache ) then
            call reaclib_interpolate_partition_functions(t9)
         end if
         call reaclib_calculate_rates(t9, rho, ye)
   end subroutine reaclib_wrapper


   subroutine reaclib_init()
         integer  :: i, ii, idx

         reaclib_ntrans(:,:,:,:) =  0
         ntarget(:,:,:,:)        =  0
         mdummy(:)               =  0
         partinter(:)            =  ONE
         nireac(:)               =  -1
         ! ^_^ masks should be initialised to .true. so that all rates are evaluated during initialisation
         reaclib_r0_mask(:)      =  .true.
         reaclib_reacv_mask(:)   =  .true.
         reaclib_partfunc_mask(:) = .true.
         rrev(:) = .false.
         is_reverse(:) = .false.
         map_reaclib_to_ppn(:)   = -1
         map_ppn_to_reaclib(:)   = -1

         select case(index_reaclib)
         case(0)
            reaclib_name   =  'BASEL'
         case(1)
            reaclib_name   =  'JINAR'
         case(2)
            reaclib_name   =  'JINAC'
         case(3)
            reaclib_name   =  'JINAV'
         end select

         if (master) then

            call reaclib_read_partition_functions()
            call reaclib_read_reaction_data()

         end if

#ifndef PPN
         call reaclib_broadcasts()
#endif

         call reaclib_create_hash()


         ! reaclib "reaction" entries are sometimes only one of several
         ! contributions to a reaction rate. Set up arrays that are
         ! one-element-per-reaction rather than one-element-per-contribution
         ! containing reactant, references and chapter, and count the total
         ! *reactions*

         ii = 0
         do i = 1, ndim_rl
            if (mdummy(i) /= 1) cycle ! rate not included

            ! check for new rate (i.e. not an additional resonance contribution)
            if (.not. is_resonance(i)) then
               ii = ii + 1
               ! redefine proper parent and dougther species, proper reference.
               re1(ii)       = e1(i)
               re2(ii)       = e2(i)
               re3(ii)       = e3(i)
               re4(ii)       = e4(i)
               re5(ii)       = e5(i)
               re6(ii)       = e6(i)
               rrlreac(ii)   = rlreac(i)
               rchapter(ii)  = chapter(i)
               rrev(ii)      = (rflag2(i) == 'v')
            end if

            ! record that this contrbution, i, is part of this rate ii
            resonance_to_rate(i) = ii
         end do

         ndim_rlt = ii

         ! initialise a, z and state information for the reactants of each reaction
         ! from reaclib. Slowest part used to be calling ispe1, which gets the index
         ! of the isotope in the reaclib nuclear data arrays. I implemented a
         ! seemingly overly-complicated search algorithm in ispe1, but I could not
         ! do any faster. Now I updated it with ispe2, which uses a hash table :)

         do i = 1, ndim_rlt
            idx       = ispe2(re1(i))
            zn1(i)    = znumpart(idx)
            an1(i)    = anumpart(idx)
            ista1(i)  = istate(idx)

            idx       = ispe2(re2(i))
            zn2(i)    = znumpart(idx)
            an2(i)    = anumpart(idx)
            ista2(i)  = istate(idx)

            idx       = ispe2(re3(i))
            zn3(i)    = znumpart(idx)
            an3(i)    = anumpart(idx)
            ista3(i)  = istate(idx)

            idx       = ispe2(re4(i))
            zn4(i)    = znumpart(idx)
            an4(i)    = anumpart(idx)
            ista4(i)  = istate(idx)

            idx       = ispe2(re5(i))
            zn5(i)    = znumpart(idx)
            an5(i)    = anumpart(idx)
            ista5(i)  = istate(idx)

            idx       = ispe2(re6(i))
            zn6(i)    = znumpart(idx)
            an6(i)    = anumpart(idx)
            ista6(i)  = istate(idx)
         end do

         call reaclib_calculate_ntrans_for_merge()

   end subroutine reaclib_init



   logical function contribution_in_ppn_network(i,ii)

         ! ^_^ check whether a reaclib rate contribution should be included in
         !     the ppn rate calculations based on its reactants and the species we
         !     have in ppn's zis isotope array. If a reactant is blank, i.e.
         !     '     ', then it is considered as included unless *all* reactants
         !     are blank.
         !
         !     the choice from where the code takes the reaction rate is also checked,
         !     along with whether he reaction will even be
         !     considered at all if it is not reaclib, then we do not need to compute it, of course.

         use reaction_info, only: lab, labb
         implicit none
         integer, intent(in)  :: i,ii
         integer              :: j, reaclib_idx(6), a, z, mapmap, idx
         logical              :: is_it_in(6), is_it_blank(6)
         character(len=5)          :: ppn_name, reactants(6)

         contribution_in_ppn_network = .false.

         ! ^_^ is it even considered at all?

         idx = map_reaclib_to_ppn(ii)
         if ( idx == -1 ) return

         mapmap = locate_vbranch_in_reaclib_rv(idx)
         if ( mapmap == -1 ) return

         if ( .not. considerreaction(mapmap) ) return

         ! ^_^ first check whether we even take this rate from reaclib
         if ( lab(mapmap) /= reaclib_name ) return

         ! ^_^ are the reactants/products all in the network?
         is_it_in(:)    = .false.
         is_it_blank(:) = .false.
         reactants(:) = [e1(i), e2(i), e3(i), e4(i), e5(i), e6(i)]

         do j = 1, 6
            if (reactants(j) == '     ') then
               is_it_in(j)    = .true.
               is_it_blank(j) = .true.
               cycle
            end if
            reaclib_idx(j) = ispe2(reactants(j))
            a = anumpart(reaclib_idx(j))
            z = znumpart(reaclib_idx(j))
            ppn_name = epsi(a,z)
            if (ppn_name /= 'nojoy') then
               is_it_in(j) = .true.
            end if
         end do

         ! this can be an issue where we assume instantaneous decay at the edges of the
         ! network in natashamaclone, so commented out for now; instead we count the
         ! contribution as in if *any* of the reactants are in the ppn network
         !if (all(is_it_in)) then
         if (any(is_it_in)) then
            if (.not. all(is_it_blank)) then
               contribution_in_ppn_network = .true.
            end if
         end if

   end function contribution_in_ppn_network


   !> create masking arrays that indicate the reactions/contributions
   !> from reaclib that we actually want to calculate for the ppn
   !> network.
   subroutine reaclib_create_masks()
#ifndef PPN
         use communication
#endif
         integer :: i, ii, a, z

         if (master) print *, " create reaclib masks ... "

         reaclib_r0_mask(:)       = .false.
         reaclib_reacv_mask(:)    = .false.
         reaclib_partfunc_mask(:) = .false.

         ii = 0

         do i = 1, ndim_rl
            if (mdummy(i) /= 1) cycle ! rate not included

            ii = resonance_to_rate(i)
            if (contribution_in_ppn_network(i,ii)) then
               reaclib_r0_mask(i)      = .true.
               reaclib_reacv_mask(ii)  = .true.
            end if
         end do

         do i = 1, nnpartdim
            a = anumpart(i)
            z = znumpart(i)
            if ( species_in_network(a,z) ) reaclib_partfunc_mask(i) = .true.
         end do
   end subroutine reaclib_create_masks


   function species_in_network( a, z ) 
         integer :: a, z
         character(len=5) :: ppn_name
         logical species_in_network

         ppn_name = epsi(a,z)
         if (ppn_name /= 'nojoy') then
            species_in_network = .true.
         else
            species_in_network = .false.
         end if
   end function species_in_network


   !> N_A<sv> (or rate(s^-1) for beta decays and photodisintegrations)
   !> if there are resonance contributions, they are calculated and
   !> are added to the main rate.
   subroutine reaclib_calculate_rates(t9in, rhoin, yein)
         use utils, only: sortedsum
         use reaction_info
         real(r8), intent(in) :: rhoin, yein, t9in
         real(r8) :: t9, t9inv, t9pow13, t9pow13inv, t9pow53, t9pow23, t9pow43, lnt9, rhoye, &
               Gfac, g1, g2, g3, g4, g5, g6, &
               t9_lower, Dt9, terms(7)
         integer  :: i, j, ii, i1, i2, i3, i4, i5, i6
         ! debugging:
         integer :: debug_fh


         t9          = max(t9in, t9_lowlimit_reaclib)
         rhoye       = rhoin * yein

         if ( use_cache .and. cache_exists ) then
            ! interpolate the pre-processed rates
            i1 = floor( ( t9 - t9_lowlimit_reaclib ) / dt9_preproc ) + 1
            i1 = max( min( i1, num_t9_preproc - 1 ), 1 )
            i2 = i1 + 1 ; t9_lower = t9_preproc(i1)
            Dt9 = ( t9 - t9_lower ) / dt9_preproc
            drates(:) = reaclib_cache(:,i2) - reaclib_cache(:,i1)
            do j = 1, jtotdim
               i = map_ppn_to_reaclib(j)
               if (.not. reaclib_reacv_mask(i)) cycle

               rv(j) = reaclib_cache(j,i1) + Dt9 * drates(j)

               ! ^_^ account for ye and density dependencies
               select case(nireac(j))
               case(i_ec)
                  ! b+, ec
                  if ( rrlreac(i) == '  ec' ) then
                     rv(j) = rv(j) * rhoye
                  end if
               case(i_ng,i_pg,i_ag,i_np,i_na,i_pn,i_pa,i_an,i_ap)
                  ! (n,g), (p,g), (a,g), (n,p), (n,a), (p,n), (p,a), (a,n), (a,p)
                  rv(j) = rv(j) * rhoin
               end select
            end do
         else
            ! ^_^ calculate rates
            t9inv       = ONE / t9
            t9pow13     = t9 ** THIRD
            t9pow13inv  = ONE / t9pow13
            t9pow23     = t9pow13 * t9pow13
            t9pow43     = t9 * t9pow13
            t9pow53     = t9pow13 * t9pow43
            lnt9        = log(t9)

            reacv(:)    = ZERO
            ii          = 0

            do i = 1, ndim_rl
               if ( .not. reaclib_r0_mask(i) ) cycle
               if ( mdummy(i) /= 1 ) cycle ! rate not included
               r0(i) = &
                     a0(i) + &
                     a1(i) * t9inv + &
                     a2(i) * t9pow13inv + &
                     a3(i) * t9pow13 + &
                     a4(i) * t9 + &
                     a5(i) * t9pow53 + &
                     a6(i) * lnt9

               r0(i) = exp(r0(i))

               ii = resonance_to_rate(i)
               reacv(ii) = reacv(ii) + r0(i)
            end do

            rv = 1.e-99_r8

            do j = 1, jtotdim

               i = map_ppn_to_reaclib(j)

               if (.not. reaclib_reacv_mask(i)) cycle

               Gfac = ONE

               if (rflag2(i) == 'v') then
                  ! calculate partition function factors and multiply the reverse rates by them.

                  i1 = ipart(an1(i),zn1(i),ista1(i))
                  i2 = ipart(an2(i),zn2(i),ista2(i))
                  i3 = ipart(an3(i),zn3(i),ista3(i))
                  i4 = ipart(an4(i),zn4(i),ista4(i))
                  i5 = ipart(an5(i),zn5(i),ista5(i))
                  i6 = ipart(an6(i),zn6(i),ista6(i))

                  g1 = partinter(i1)
                  g2 = partinter(i2)
                  g3 = partinter(i3)
                  g4 = partinter(i4)
                  g5 = partinter(i5)
                  g6 = partinter(i6)

                  select case(rchapter(i))
                  case('1')
                     Gfac = g2 / g1
                  case('2')
                     Gfac = g2 * g3 / g1
                  case('3')
                     Gfac = g2 * g3 * g4 / g1
                  case('4')
                     Gfac = g3 / (g1 * g2)
                  case('5')
                     Gfac = g3 * g4 / (g1 * g2)
                  case('6')
                     Gfac = g3 * g4 * g5 / (g1 * g2)
                  case('7')
                     Gfac = g3 * g4 * g5 * g6 / (g1 * g2)
                  case('8')
                     Gfac = g4 / (g1 * g2 * g3)
                  case('9')
                     Gfac = g4 * g5 / (g1 * g2 * g3)
                  case('10')
                     Gfac = g5 * g6 / (g1 * g2 * g3 * g4)
                  case('11')
                     Gfac = g2 * g3 * g4 * g5 * g6 / g1
                  case default
                     Gfac = ONE
                  end select
               end if

               reacv(i) = reacv(i) * Gfac

               if ( preprocessing ) then
                  rv(j) = reacv(i)
               else
                  select case(nireac(j))
                     ! RJS 23/11/19 -- treat i_bp cases like i_gp
                  case(i_bm,i_gn,i_bn,i_bp,i_gp,i_ga)
                     rv(j) = reacv(i)
                  case(i_ec)
                     if ( rrlreac(i) == '  ec' ) then
                        rv(j) = reacv(i) * rhoye
                     else
                        rv(j) = reacv(i)
                     end if
                  case(i_ng,i_pg,i_ag,i_np,i_na,i_pn,i_pa,i_an,i_ap)
                     rv(j) = reacv(i) * rhoin
                  case default
                     rv(j) = 1.e-99_r8
                  end select
               end if
            end do
         end if

         rv = max(rv, 1e-99_r8)
   end subroutine reaclib_calculate_rates


   !> Interpolate reaclib partition functions for the given temperature
   subroutine reaclib_interpolate_partition_functions(t9_in)
         use utils, only: bsearch_r8
         implicit none
         real(r8), intent(in) :: t9_in
         real(r8) :: t9, t90, t91, Dt9
         integer  :: i, i0, i1, k

         ! ^_^ clip temperature
         t9 = max(rt9(1), min(t9_in, rt9(partdim)))

         ! ^_^ bisection search for temperature
         call bsearch_r8(rt9, t9, i0)
         i1 = i0 + 1; t90 = rt9(i0); t91 = rt9(i1)
         Dt9 = ( t9 - t90 ) / ( t91 - t90 )

         do i = 1, nnpartdim
            if ( .not. reaclib_partfunc_mask(i) ) cycle
            partinter    (i) = partnum(i,i0) + (partnum(i,i1) - partnum(i,i0)) * Dt9
         end do
   end subroutine reaclib_interpolate_partition_functions



   subroutine reaclib_preprocessor
         ! ^_^ evaluate reaclib formulae for reactions that will be used in the current
         ! computation over num_t9_preproc temperature points, so that interpolation can be used
         ! during run time instead of formula evaluation
         real(r8) :: t9, tti1
         integer :: i

         preprocessing = .true.
         allocate( reaclib_cache( jtotdim, num_t9_preproc ) )
         allocate( drates(jtotdim) )
         dt9_preproc = ( 10._r8 - t9_lowlimit_reaclib ) / ( num_t9_preproc - 1 )

#ifndef PPN
         if ( master ) print *, ' creating reaclib cache...'
#else
         print *, ' creating reaclib cache...'
#endif
         tti1 = wallclocktime()

         do i = 1, num_t9_preproc
            t9 = t9_lowlimit_reaclib + ( i - 1 ) * dt9_preproc
            t9_preproc(i) = t9
            call reaclib_interpolate_partition_functions( t9 )
            call reaclib_calculate_rates( t9, ONE, ONE )
            reaclib_cache( :, i ) = rv(1:jtotdim)
            write(*,"(a17,2(es23.15))") 'reaclib preproc: ', t9, reaclib_cache( 57105, i )
         end do

         cache_exists = .true.
         preprocessing = .false.

#ifndef PPN
         if ( master ) write(*,fmtprofile) 'reaclib cache time / s = ', wallclocktime() - tti1
#else
         write(*,fmtprofile) 'reaclib cache time / s = ', wallclocktime() - tti1
#endif

   end subroutine reaclib_preprocessor



   subroutine reaclib_calculate_ntrans_for_merge()
         use reaction_info
         ! ^_^ this subroutine populates the 4d array reaclib_ntrans that is used
         !     to merge the reaclib rates into the ppn network.

         ! ^_^ previously this was called every physics call because there are some
         !     multiplications of the rate by rho or rhoye. This has been moved to the
         !     calculate_rates subroutine

         integer :: i, j

         j = 0

         do i = 1, ndim_rlt
            if (re5(i) /= '     ') cycle
            if (rchapter(i) == '1') then
               if (zn1(i) == zn2(i)-1) then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista1(i)
                  nireac(j)   = i_bm       !beta-
                  aan(j)      = an1(i)
                  zzn(j)      = zn1(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista2(i)
                  aanf(j)     = an2(i)
                  zznf(j)     = zn2(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (zn1(i) == zn2(i)+1) then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista1(i)
                  nireac(j)   = i_ec       !beta+,ec
                  aan(j)      = an1(i)
                  zzn(j)      = zn1(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista2(i)
                  aanf(j)     = an2(i)
                  zznf(j)     = zn2(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               end if
            else if (rchapter(i) == '2') then
               ! RJS 23/11/19 -- Why are we checking that the 4th slot in the reaclib list is blank?
               ! It ought to be BY DEFINITION of a reaclib chapter 2 rate, e_1 -> e_2 + e_3
               if (re2(i) == '    n' .and. re4(i) == '     ') then
                  if (zn1(i) == zn3(i)) then
                     j = j + 1
                     map_reaclib_to_ppn(i) = j
                     map_ppn_to_reaclib(j) = i
                     ista(j)    = ista1(i)
                     nireac(j)  = i_gn       !g,n
                     aan(j)     = an1(i)
                     zzn(j)     = zn1(i)
                     reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                     istaf(j)   = ista3(i)
                     aanf(j)    = an3(i)
                     zznf(j)    = zn3(i)
                     ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                     is_reverse(j) = rrev(i)
                  else if (zn1(i)+1 == zn3(i)) then
                     j = j + 1
                     map_reaclib_to_ppn(i) = j
                     map_ppn_to_reaclib(j) = i
                     ista(j)    = ista1(i)
                     nireac(j)  = i_bn       !beta-delayed neutron emission
                     aan(j)     = an1(i)
                     zzn(j)     = zn1(i)
                     reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                     istaf(j)   = ista3(i)
                     aanf(j)    = an3(i)
                     zznf(j)    = zn3(i)
                     ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                     is_reverse(j) = rrev(i)
                  end if
               else if (re2(i) == '    p' .and. re4(i) == '     ') then
                  ! RJS 23/11/19 -- adding test to distinguish between g,p and b,p reactions
                  ! Use test based on proton number
                  if (zn1(i) - 1 == zn3(i) ) then
                     j = j + 1
                     map_reaclib_to_ppn(i) = j
                     map_ppn_to_reaclib(j) = i
                     ista(j)     = ista1(i)
                     nireac(j)   = i_gp       !g,p
                     aan(j)      = an1(i)
                     zzn(j)      = zn1(i)
                     reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                     istaf(j)    = ista3(i)
                     aanf(j)     = an3(i)
                     zznf(j)     = zn3(i)
                     ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                     is_reverse(j) = rrev(i)
                  else if (zn1(i) - 2 == zn3(i)) then
                      j = j + 1
                      map_reaclib_to_ppn(i) = j
                      map_ppn_to_reaclib(j) = i
                      ista(j)     = ista1(i)
                      nireac(j) = i_bp   ! beta-delayed proton emission
                      aan(j)      = an1(i)
                      zzn(j)      = zn1(i)
                      reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                      istaf(j)    = ista3(i)
                      aanf(j)     = an3(i)
                      zznf(j)     = zn3(i)
                      ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                      is_reverse(j) = rrev(i)
                   end if
               else if (re2(i) == '  he4' .and. re4(i) == '     ') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista1(i)
                  nireac(j)   = i_ga      !g,a
                  aan(j)      = an1(i)
                  zzn(j)      = zn1(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista3(i)
                  aanf(j)     = an3(i)
                  zznf(j)     = zn3(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               end if
            else if (rchapter(i) == '4') then
               if (re1(i) == '    n' .and. re4(i) == '     ') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_ng       !n,g
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista3(i)
                  aanf(j)     = an3(i)
                  zznf(j)     = zn3(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '    p' .and. re4(i) == '     ') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_pg       !p,g
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista3(i)
                  aanf(j)     = an3(i)
                  zznf(j)     = zn3(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '  he4' .and. re4(i) == '     ') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_ag       !a,g
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista3(i)
                  aanf(j)     = an3(i)
                  zznf(j)     = zn3(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               end if
            else if (rchapter(i) == '5') then
               if (re1(i) == '    n' .and. re3(i) == '    p') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_np       !n,p
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '    n' .and. re3(i) == '  he4') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_na       !n,a
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '    p' .and. re3(i) == '    n') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_pn       !p,n
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '    p' .and. re3(i) == '  he4') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_pa       !p,a
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '  he4' .and. re3(i) == '    n') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_an      !a,n
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               else if (re1(i) == '  he4' .and. re3(i) == '    p') then
                  j = j + 1
                  map_reaclib_to_ppn(i) = j
                  map_ppn_to_reaclib(j) = i
                  ista(j)     = ista2(i)
                  nireac(j)   = i_ap      !a,p
                  aan(j)      = an2(i)
                  zzn(j)      = zn2(i)
                  reaclib_ntrans(aan(j),zzn(j),ista(j),nireac(j)) = j
                  istaf(j)    = ista4(i)
                  aanf(j)     = an4(i)
                  zznf(j)     = zn4(i)
                  ntarget(aanf(j),zznf(j),istaf(j),nireac(j)) = j
                  is_reverse(j) = rrev(i)
               end if
            end if
         end do

         jtotdim = j

   end subroutine reaclib_calculate_ntrans_for_merge


   subroutine reaclib_read_reaction_data()
         character(len=128)  :: reacfile
         character(len=1)    :: rcdummy
         integer        :: i, j

         select case(index_reaclib)
         case(0)
            reacfile = '../NPDATA/REACLIB/reaclib.nosmo'
         case(1)
            reacfile = '../NPDATA/REACLIB/20081109ReaclibV0.5'
         case(2)
            !reacfile = '../NPDATA/REACLIB/20120510ReaclibV1.1'
            reacfile = '../NPDATA/REACLIB/results07010946'
         case(3)
            reacfile = '../NPDATA/REACLIB/results01111258'
         end select

         open(unit=table_input, file = reacfile, status='old')

         i = 1
         j = 1
         rcdummy = ' '

         13 read(table_input, '(A1)', end=992) cdummy

         if (cdummy.eq.'#') then
            goto 992
         else if (cdummy.eq.' ') then
            backspace(table_input)
            ! *** command to skip elements beyond At in reaclib jina.
            ! *** not avaliable for now max excess, partition function, etc.
            read(table_input, '(5x,5A5)', end=992) cdummy1,cdummy2,cdummy3,cdummy4,cdummy5
            if (out_of_bounds()) then
               read(table_input,*)
               read(table_input,*)
               goto 13
            else
               backspace(table_input)
            end if
            read(table_input, '(A1,4x,6A5,8x,A4,A1,A1,3x,1ES12.5)', end=992) &
                  cdummy, e1(i), e2(i), e3(i), e4(i), e5(i), e6(i), &
                  rlreac(i), rflag1(i), rflag2(i), qmev(i)
            read(table_input,'(4ES13.6)') a0(i), a1(i), a2(i), a3(i)
            read(table_input,'(3ES13.6)') a4(i), a5(i), a6(i)
            chapter(i) = rcdummy
            mdummy(i) = 1
            i = i + 1
            goto 13
         else
            chapter(i)  = cdummy
            rcdummy     = cdummy
            read(table_input,*)
            read(table_input,*)
            j = j + 1
            goto 13
         endif

         992 continue

         close(table_input)

   end subroutine reaclib_read_reaction_data



   subroutine reaclib_read_partition_functions()
         use physics_knobs, only: index_reaclib
         implicit none
         real(r8) :: aanumpart(nnpartdim)
         character(len=5)    :: ccdummy
         character(len=128)  :: partfile
         integer        :: i, k, jj

         select case(index_reaclib)
         case(0)
            partfile = '../NPDATA/REACLIB/winvn.nosmo.correct'
         case(1)
            partfile = '../NPDATA/REACLIB/winvn.nosmo.correct_jina'
         case(2)
            !partfile = '../NPDATA/REACLIB/winvne_v1.0_modified.dat'
            partfile = '../NPDATA/REACLIB/winvne_v9.9_modified.dat'
         case(3)
            partfile = '../NPDATA/REACLIB/winvne_v1.0_modified_cf.dat'
         end select

         open(unit=table_input, file = partfile, status='old')

         do i = 1, 2
            read(table_input,*)
         end do

         read(table_input,*) (rt9(k), k = 1, partdim)
         11 read(table_input,'(4x,a1)',end=990) cdummy
         if (cdummy.ne.'n') then
            goto 11
         else
            backspace(table_input)
         endif
         990     continue

         i = 1
         12 read(table_input,'(a1)',end=991)ccdummy
         if (ccdummy.ne.'#') then
            backspace(table_input)
            ccdummy = spe_part(i)
            read(table_input,'(A5,F12.3,2(I4),F6.1,F10.3)', end=991) &
                  spe_part(i), aanumpart(i), znumpart(i), nnumpart(i), &
                  spinpart(i), exmass(i)
            read(table_input,*) (partnum(i,jj), jj=1,partdim)
            anumpart(i) = int(aanumpart(i))

            ! check if there is ground state + isomeric state or just thermalized
            ! state.
            ! define index ipart to identify the species.

            if (spe_part(i)(4:4) == '-') then
               istate(i) = 2
            else if (spe_part(i)(4:4) == '*') then
               istate(i) = 3
            else
               istate(i) = 1
            end if

            ipart(anumpart(i),znumpart(i),istate(i)) = i

            i = i + 1

            goto 12
         else
            goto 991
         endif

         991 continue
         close(table_input)

         ! *** definition of ipart for empty specie and for its partition function.
         i = nnpartdim + 1
         select case(index_reaclib)
         case(3)
            anumpart(i)    = i325dim+1
            znumpart(i)    = iCfdim+1
         case default
            anumpart(i)    = i282dim+1
            znumpart(i)    = iAtdim+1
         end select
         istate(i)      = 1
         partinter(i)   = 1._r8
         spe_part(i)    = '     '
         ipart(anumpart(i),znumpart(i),istate(i)) = i

   end subroutine reaclib_read_partition_functions



   integer function ispe2(name1)

         ! ^_^ an upgrade of ispe1, the function that returns the index of a reaclib
         !     isotope name in the reaclib nuclear data arrays. This version does a
         !     lookup of a hash table that was created during reaclib_init

         character(len=5), intent(in) :: name1
         integer                 :: idx, i

         if (.not. allocated(reaclib_name_to_idx% values)) then
            stop 'attempting to access reaclib hash array before it has been created'
         end if

         ispe2 = -1
         idx = hash(name1)
         do i = 1, reaclib_name_to_idx% nvalsinbucket
            if (reaclib_name_to_idx% values(idx,i) == -1) then
               print *, 'error: isotope not found in reaclib hash table', name1
               stop
            else if (name1 == spe_part(reaclib_name_to_idx% values(idx,i))) then
               ispe2 = reaclib_name_to_idx% values(idx,i)
               exit
            end if
         end do

   end function ispe2


   integer function ispe1(name1)

         integer, parameter :: nspe = 101
         integer ifoundnumber,i,A,inum,inam,idx,cou,cnt,Z, DZ(nspe), idx1, breakout
         character(len=5) name1,nam,num,d(nspe)
         character(len=1) ch

         ispe1 = -1 ! assume this value if no valid name is found

         ifoundnumber = 0

         D = [ 'n0', 'h ', 'd ', 't ', 'he', 'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', &
               'ne', 'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', 'k ', 'ca', &
               'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', &
               'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y ', 'zr', 'nb', 'mo', &
               'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i ', &
               'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', &
               'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', &
               'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', &
               'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu', 'am', 'cm', 'bk', &
               'cf' ]

         DZ = [ 0,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, &
               10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
               21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
               32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, &
               43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, &
               54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, &
               65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, &
               76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
               87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
               98 ]

         ! *** with the if statements below I am skipping the do loop in case the specie
         ! *** is a neutron, a proton, he4 (position 1,2,6 in the species list--> no
         !-->time saved).
         if(name1 .eq. '    n') then
            ispe1 = 1
            ifoundnumber = 1
            goto 11
         else if(name1 .eq. '    p') then
            ispe1 = 2
            ifoundnumber = 1
            goto 11
         else if(name1 .eq. '  he4') then
            ispe1 = 6
            ifoundnumber = 1
            goto 11
         else if(name1 .eq. '     ') then
            ispe1 = nnpartdim + 1
            ifoundnumber = 1
            goto 11
         end if

         !! seperate name and number
         inum=1
         inam=1
         do i=1,5
            nam(i:i) = ' '
            num(i:i) = ' '
         end do
         do i=1,5
            ch = NAME1(i:i)
            if ( is_numeric(ch) ) then
               num(inum:inum) = ch
               inum = inum + 1
            else if ( ch .eq. ' ') then
               cycle
            else
               nam(inam:inam) = ch
               inam = inam + 1
            end if
         end do

         Z=-1
         do i = 1, nspe
            if ( trim(nam) .eq. trim(D(i)) ) then
               Z = DZ(i)
               exit
            end if
         end do

         !! search spe_part for the index
         idx = Z * ( Z + 2 ) ! sqrt(real(Z)) * ( Z + 2 )
         idx = max(min( idx, nnpartdim+1 ),1)
         idx1=idx
         cou = 1
         cnt = 1
         if ( trim(nam) .eq. 'al-' .or. trim(nam) .eq. 'al*' ) then
            do I=1,nnpartdim+1
               if(name1 .eq. spe_part(i)) then
                  ispe1=i
                  ifoundnumber = 1
                  exit
               endif
            end do
         else
            !!----------------
            !! this one starts with the educated guess and then looks
            !! left, right, left, right etc to find the index.
            do
               if ( cnt .gt. nnpartdim+1 ) exit
               if ( spe_part(idx1) .eq. name1 ) then
                  ispe1 = idx1
                  ifoundnumber = 1
                  exit
               else
                  idx = idx + (-1) ** mod(cnt,2) * cou
                  idx1 = min(max(idx,1), nnpartdim+1)
                  cou = cou + 1
               end if
               cnt = cnt + 1
            end do
         end if

         11    continue

   end function ispe1



   logical function is_numeric(string)
         character(len=*), intent(in) :: string
         real :: x
         integer :: e

         read(string,*,iostat=e) x
         is_numeric = e == 0
   end function is_numeric



#ifndef PPN
   subroutine reaclib_broadcasts()

         ! ^_^ integers
         call broadcast(mdummy)
         call broadcast(ipart)
         call broadcast(istate)
         call broadcast(anumpart)
         call broadcast(znumpart)

         ! ^_^ r8 reals
         call broadcast(rt9)
         call broadcast(exmass)
         call broadcast(qmev)
         call broadcast(a0)
         call broadcast(a1)
         call broadcast(a2)
         call broadcast(a3)
         call broadcast(a4)
         call broadcast(a5)
         call broadcast(a6)
         call broadcast(partnum)
         call broadcast(spinpart)

         ! ^_^ characters
         call broadcast_ch_arr(spe_part)
         call broadcast_ch_arr(rflag1)
         call broadcast_ch_arr(rflag2)
         call broadcast_ch_arr(rlreac)
         call broadcast_ch_arr(e1)
         call broadcast_ch_arr(e2)
         call broadcast_ch_arr(e3)
         call broadcast_ch_arr(e4)
         call broadcast_ch_arr(e5)
         call broadcast_ch_arr(e6)
         call broadcast_ch_arr(chapter)

   end subroutine reaclib_broadcasts
#endif
   !(PPN)


   logical function out_of_bounds()
         ! ^_^ TODO: use/create funtion returning Z given element name and then just check if Z larger than whatever...

         if ( index_reaclib < 3 ) then
            if (cdummy1(1:2).eq.'rn'.or.cdummy1(1:2).eq.'fr'.or. &
                  cdummy1(1:2).eq.'ra'.or.cdummy1(1:2).eq.'ac'.or. &
                  cdummy1(1:2).eq.'th'.or.cdummy1(1:2).eq.'pa'.or. &
                  cdummy1(1:2).eq.' u'.or.cdummy1(1:2).eq.'np'.or. &
                  cdummy1(1:2).eq.'pu'.or.cdummy1(1:2).eq.'am'.or. &
                  cdummy1(1:2).eq.'cm'.or.cdummy1(1:2).eq.'bk'.or. &
                  cdummy1(1:2).eq.'cf'.or.cdummy1(1:2).eq.'es'.or. &
                  cdummy1(1:2).eq.'fm'.or.cdummy1(1:2).eq.'md'.or. &
                  cdummy1(1:2).eq.'no'.or.cdummy1(1:2).eq.'lr'.or. &
                  cdummy1(1:2).eq.'rf'.or.cdummy1(1:2).eq.'db'.or. &
                  cdummy1(1:2).eq.'sg'.or.cdummy1(1:2).eq.'bh'.or. &
                  cdummy1(1:2).eq.'hs'.or.cdummy1(1:2).eq.'mt'.or. &
                  cdummy1(1:2).eq.'ds'.or.cdummy1(1:2).eq.'rg'.or. &
                  cdummy2(1:2).eq.'rn'.or.cdummy2(1:2).eq.'fr'.or. &
                  cdummy2(1:2).eq.'ra'.or.cdummy2(1:2).eq.'ac'.or. &
                  cdummy2(1:2).eq.'th'.or.cdummy2(1:2).eq.'pa'.or. &
                  cdummy2(1:2).eq.' u'.or.cdummy2(1:2).eq.'np'.or. &
                  cdummy2(1:2).eq.'pu'.or.cdummy2(1:2).eq.'am'.or. &
                  cdummy2(1:2).eq.'cm'.or.cdummy2(1:2).eq.'bk'.or. &
                  cdummy2(1:2).eq.'cf'.or.cdummy2(1:2).eq.'es'.or. &
                  cdummy2(1:2).eq.'fm'.or.cdummy2(1:2).eq.'md'.or. &
                  cdummy2(1:2).eq.'no'.or.cdummy2(1:2).eq.'lr'.or. &
                  cdummy2(1:2).eq.'rf'.or.cdummy2(1:2).eq.'db'.or. &
                  cdummy2(1:2).eq.'sg'.or.cdummy2(1:2).eq.'bh'.or. &
                  cdummy2(1:2).eq.'hs'.or.cdummy2(1:2).eq.'mt'.or. &
                  cdummy2(1:2).eq.'ds'.or.cdummy2(1:2).eq.'rg'.or. &
                  cdummy3(1:2).eq.'rn'.or.cdummy3(1:2).eq.'fr'.or. &
                  cdummy3(1:2).eq.'ra'.or.cdummy3(1:2).eq.'ac'.or. &
                  cdummy3(1:2).eq.'th'.or.cdummy3(1:2).eq.'pa'.or. &
                  cdummy3(1:2).eq.' u'.or.cdummy3(1:2).eq.'np'.or. &
                  cdummy3(1:2).eq.'pu'.or.cdummy3(1:2).eq.'am'.or. &
                  cdummy3(1:2).eq.'cm'.or.cdummy3(1:2).eq.'bk'.or. &
                  cdummy3(1:2).eq.'cf'.or.cdummy3(1:2).eq.'es'.or. &
                  cdummy3(1:2).eq.'fm'.or.cdummy3(1:2).eq.'md'.or. &
                  cdummy3(1:2).eq.'no'.or.cdummy3(1:2).eq.'lr'.or. &
                  cdummy3(1:2).eq.'rf'.or.cdummy3(1:2).eq.'db'.or. &
                  cdummy3(1:2).eq.'sg'.or.cdummy3(1:2).eq.'bh'.or. &
                  cdummy3(1:2).eq.'hs'.or.cdummy3(1:2).eq.'mt'.or. &
                  cdummy3(1:2).eq.'ds'.or.cdummy3(1:2).eq.'rg'.or. &
                  cdummy4(1:2).eq.'rn'.or.cdummy4(1:2).eq.'fr'.or. &
                  cdummy4(1:2).eq.'ra'.or.cdummy4(1:2).eq.'ac'.or. &
                  cdummy4(1:2).eq.'th'.or.cdummy4(1:2).eq.'pa'.or. &
                  cdummy4(1:2).eq.' u'.or.cdummy4(1:2).eq.'np'.or. &
                  cdummy4(1:2).eq.'pu'.or.cdummy4(1:2).eq.'am'.or. &
                  cdummy4(1:2).eq.'cm'.or.cdummy4(1:2).eq.'bk'.or. &
                  cdummy4(1:2).eq.'cf'.or.cdummy4(1:2).eq.'es'.or. &
                  cdummy4(1:2).eq.'fm'.or.cdummy4(1:2).eq.'md'.or. &
                  cdummy4(1:2).eq.'no'.or.cdummy4(1:2).eq.'lr'.or. &
                  cdummy4(1:2).eq.'rf'.or.cdummy4(1:2).eq.'db'.or. &
                  cdummy4(1:2).eq.'sg'.or.cdummy4(1:2).eq.'bh'.or. &
                  cdummy4(1:2).eq.'hs'.or.cdummy4(1:2).eq.'mt'.or. &
                  cdummy4(1:2).eq.'ds'.or.cdummy4(1:2).eq.'rg'.or. &
                  cdummy5(1:2).eq.'rn'.or.cdummy5(1:2).eq.'fr'.or. &
                  cdummy5(1:2).eq.'ra'.or.cdummy5(1:2).eq.'ac'.or. &
                  cdummy5(1:2).eq.'th'.or.cdummy5(1:2).eq.'pa'.or. &
                  cdummy5(1:2).eq.' u'.or.cdummy5(1:2).eq.'np'.or. &
                  cdummy5(1:2).eq.'pu'.or.cdummy5(1:2).eq.'am'.or. &
                  cdummy5(1:2).eq.'cm'.or.cdummy5(1:2).eq.'bk'.or. &
                  cdummy5(1:2).eq.'cf'.or.cdummy5(1:2).eq.'es'.or. &
                  cdummy5(1:2).eq.'fm'.or.cdummy5(1:2).eq.'md'.or. &
                  cdummy5(1:2).eq.'no'.or.cdummy5(1:2).eq.'lr'.or. &
                  cdummy5(1:2).eq.'rf'.or.cdummy5(1:2).eq.'db'.or. &
                  cdummy5(1:2).eq.'sg'.or.cdummy5(1:2).eq.'bh'.or. &
                  cdummy5(1:2).eq.'hs'.or.cdummy5(1:2).eq.'mt'.or. &
                  cdummy5(1:2).eq.'ds'.or.cdummy5(1:2).eq.'rg') then
            out_of_bounds = .true.
         else
            out_of_bounds = .false.
         end if
      else ! index_reaclib == 3
         if (cdummy1(1:2).eq.'es'.or. &
               cdummy1(1:2).eq.'fm'.or.cdummy1(1:2).eq.'md'.or. &
               cdummy1(1:2).eq.'no'.or.cdummy1(1:2).eq.'lr'.or. &
               cdummy1(1:2).eq.'rf'.or.cdummy1(1:2).eq.'db'.or. &
               cdummy1(1:2).eq.'sg'.or.cdummy1(1:2).eq.'bh'.or. &
               cdummy1(1:2).eq.'hs'.or.cdummy1(1:2).eq.'mt'.or. &
               cdummy1(1:2).eq.'ds'.or.cdummy1(1:2).eq.'rg'.or. &
               cdummy2(1:2).eq.'es'.or. &
               cdummy2(1:2).eq.'fm'.or.cdummy2(1:2).eq.'md'.or. &
               cdummy2(1:2).eq.'no'.or.cdummy2(1:2).eq.'lr'.or. &
               cdummy2(1:2).eq.'rf'.or.cdummy2(1:2).eq.'db'.or. &
               cdummy2(1:2).eq.'sg'.or.cdummy2(1:2).eq.'bh'.or. &
               cdummy2(1:2).eq.'hs'.or.cdummy2(1:2).eq.'mt'.or. &
               cdummy2(1:2).eq.'ds'.or.cdummy2(1:2).eq.'rg'.or. &
               cdummy3(1:2).eq.'es'.or. &
               cdummy3(1:2).eq.'fm'.or.cdummy3(1:2).eq.'md'.or. &
               cdummy3(1:2).eq.'no'.or.cdummy3(1:2).eq.'lr'.or. &
               cdummy3(1:2).eq.'rf'.or.cdummy3(1:2).eq.'db'.or. &
               cdummy3(1:2).eq.'sg'.or.cdummy3(1:2).eq.'bh'.or. &
               cdummy3(1:2).eq.'hs'.or.cdummy3(1:2).eq.'mt'.or. &
               cdummy3(1:2).eq.'ds'.or.cdummy3(1:2).eq.'rg'.or. &
               cdummy4(1:2).eq.'es'.or. &
               cdummy4(1:2).eq.'fm'.or.cdummy4(1:2).eq.'md'.or. &
               cdummy4(1:2).eq.'no'.or.cdummy4(1:2).eq.'lr'.or. &
               cdummy4(1:2).eq.'rf'.or.cdummy4(1:2).eq.'db'.or. &
               cdummy4(1:2).eq.'sg'.or.cdummy4(1:2).eq.'bh'.or. &
               cdummy4(1:2).eq.'hs'.or.cdummy4(1:2).eq.'mt'.or. &
               cdummy4(1:2).eq.'ds'.or.cdummy4(1:2).eq.'rg'.or. &
               cdummy5(1:2).eq.'es'.or. &
               cdummy5(1:2).eq.'fm'.or.cdummy5(1:2).eq.'md'.or. &
               cdummy5(1:2).eq.'no'.or.cdummy5(1:2).eq.'lr'.or. &
               cdummy5(1:2).eq.'rf'.or.cdummy5(1:2).eq.'db'.or. &
               cdummy5(1:2).eq.'sg'.or.cdummy5(1:2).eq.'bh'.or. &
               cdummy5(1:2).eq.'hs'.or.cdummy5(1:2).eq.'mt'.or. &
               cdummy5(1:2).eq.'ds'.or.cdummy5(1:2).eq.'rg') then
         out_of_bounds = .true.
      else
         out_of_bounds = .false.
      end if
   end if

   end function out_of_bounds



   logical function is_resonance(i)

         ! ^_^ is the reaclib entry at i a resonant contribution to the previous entry
         !     or is it a new reaction entirely?

         integer, intent(in) :: i


         if (i == 1) then
            is_resonance = .false.
         else if (e1(i) /= e1(i-1) .or. e2(i) /= e2(i-1) .or. &
               e3(i) /= e3(i-1) .or. e4(i) /= e4(i-1) .or. &
               e5(i) /= e5(i-1) .or. e6(i) /= e6(i-1)) then
            is_resonance = .false.
         else
            is_resonance = .true.
         end if

   end function is_resonance



end module reaclib

