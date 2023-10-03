! define PROFILE to print profiling information
#undef PROFILE

module physics
   use communication
   use utils, only: r8, i4
   implicit none
   private

   public rnetw2007, rnetw2008

contains

   !> This subroutine creates the network, including defining all isotopes and isomers and the
   !> reactions between them.
   !> 
   !> output:
   !> vbranch(nsp,num_rtypes) : reaction rates for this channel
   !> nbranch(nsp,num_rtypes) : target (product) of reaction, adjusted for network boundaries
   !> label  (nsp,num_rtypes) : reference for the reaction rate
   !
   !  NOTE from Marco:
   !  amin(isotot) is a parameter used in case of network
   !  uncomplete.  if a reaction produces an unstable isotope not
   !  included in the network, an istantaneous decay is assumed. beta+
   !  or beta-? anim(isotot) helps in this way, because for every
   !  element is fixed an amin(isotot). if a(unstable isotope) <
   !  amin(element) then  beta+ decay is assumed, if a(unstable isotope)
   !  > amin(element) then  beta- decay is assumed. in most of the cases
   !  amin(el) is the lighter stable isotope of element. few exceptions
   !  are present: EU152,TB158,TA180,AU196 (because these isotopes are
   !  ustable for both beta- and beta+) and tc98 and PM146 (no stable
   !  isotopes are present for these elements).  nb: before using a new
   !  network this part should be controlled again in order to avoid
   !  mistakes!!!!!!  nb: beta+ unstable isotopes heavier than the
   !  lightest stable isotope should be included!!
   !
   !  NB: starting from version 8 amin(el) is not used! it is used a kind of
   !  y_e for beta+/e+ region.
   subroutine networkI ( ye, t9, rho, nbranch, vbranch, speneu, label )
         use array_sizes
         use reaction_info, only: num_rtypes, i_ng, i_gn, i_np, i_na, i_pg, i_gp, i_pn, i_pa, &
              ! RJS 23/11/19 -- added i_bp
               i_ag, i_ga, i_an, i_ap, i_bm, i_ec, i_bn, i_bp
         use netgen, only: nacre_name, netgen_beta_ntrans, netgen_nacre_ntrans, &
               netgen_illi_ntrans, netgen_oda_ntrans, netgen_lmp_ntrans, &
               netgen_interpolate_rates, netgen_rtype, beta, illi, nacre, oda, lmp
         use fuller, only: fuller_ntrans, fuller_interpolate_rates, fuller_rates
         use physics_knobs
         use frame_knobs
         use utils
         use reaclib, only: aan, zzn, aanf, zznf, reaclib_ntrans, nireac, jtotdim, reaclib_name, &
               reaclib_wrapper, rv
         use kadonis, only: nkad, zkad, akad, vkad, kadonis_interpolate_rates
         use jbj16, only: jbj_rate, jbjminA, jbjmaxA, jbjminZ, jbjmaxZ, &
               jbj_ntrans, jbj_interpolate_rates
         use nkk04
         use network_boundaries, only: natashamcloane
         use nuc_data, only: amin, atot, ztot, ntot, eltot, read_isotopedatabase
         use other_nuc
         use constants
         use communication

         implicit none

         integer :: &
               i, j, k, kk, & !< counters for various things
               rtype, &       !< reaction type
               aareac, &      !< flag 0 or 1: b-plus or b-minus decay isotope at edge of network
               isotot, &      !< isotot is the equivalent to nnn in this routine
               nbranch(nsp,num_rtypes)

         real(r8) :: t9, ye, rho, betathresh, vbranch(nsp,num_rtypes)

         ! NETGEN
         integer :: ireac
         character(len=5) :: label(nsp,num_rtypes) !< references for the reaction rates
         character(len=5) :: speneu(nsp) !< isotope names, that will eventually become zis in nuc_data

         logical considerisoinlist(nsp)

         ! index for weak reactions full - seems to just point to the species given its a and z...
         ! so doesn't seem to be unique to weak reactions...
#if pIDX_RCLB == 3
         integer :: index_raw(0:iCfdim,i325dim)
#else
         integer :: index_raw(0:iAtdim,i282dim)
#endif

         index_raw(:,:) = 0

         isotot = nnn
         call read_isotopedatabase(isotot, considerisoinlist)

         ! get weak rates from netgen. These are used to determine the width of the beta-valley for
         ! the network setup (parameter tbetamin)
         call netgen_interpolate_rates(t9, rho, ye, beta)

         do kk = 1, isotot
            if (atot(kk) > 226 .or. ztot(kk) > 85) cycle ! max A and Z in netgen

            do rtype = i_bm, i_ec ! beta-decay and electron capture
               i = 0; i = netgen_beta_ntrans(atot(kk),ztot(kk),1,rtype)
               if (i <= 0) cycle
               vbranch(kk,rtype) = beta % rates(i)
            end do
         end do

         ! exclude unstable isotopes with half lifes < tbetamin: i.e. rate > ln2 / tbetamin
         betathresh = ln2 / max(tbetamin,1e-99_r8)

         k = 0
         do kk = 1, isotot
            if (.not. considerisoinlist(kk)) then
               isotot = isotot - 1
               cycle
            end if
            if (vbranch(kk,i_bm) > betathresh .or. vbranch(kk,i_ec) > betathresh) then
               isotot = isotot - 1
               cycle
            end if
            k = k + 1
            eltot  (k)   = eltot  (kk)
            atot   (k)   = atot   (kk)
            ztot   (k)   = ztot   (kk)
            ntot   (k)   = ntot   (kk)
            amin   (k)   = amin   (kk)

            index_raw(ztot(k),atot(k)) = k
         end do

         nnn = isotot

         eltot  (isotot+1:)    = 'FH'
         ztot   (isotot+1:)    = 0
         atot   (isotot+1:)    = 0
         ntot   (isotot+1:)    = 0
         amin   (isotot+1:)    = 0

         ! set all weak reaction rates back to zero
         vbranch(:,i_bm) = ZERO
         vbranch(:,i_ec) = ZERO

         ! get all element names in capital letters.
         do i = 1, nsp
            if (eltot(i) /= 'nn') eltot(i) = to_upper(eltot(i))
         end do

         call reaclib_wrapper(t9, rho, ye)

         ! assign rates to correct position according to the reaction type.
         do kk = 1, isotot
            do rtype = 1, num_rtypes
               ! find reaction rtype from species kk from reaclib
               i = reaclib_ntrans(atot(kk),ztot(kk),1,rtype)
               if (i <= 0) cycle

               ! write it to nasv and label arrays
               vbranch(kk,rtype) = rv(i)
               label(kk,rtype)   = reaclib_name
            end do
         end do

         ! assign targets from reaclib
         do i = 1, jtotdim
            do kk = 1, isotot
               ! search for isotope in our network
               if (ztot(kk) /= zzn(i) .or. atot(kk) /= aan(i)) cycle

               rtype = nireac(i)
               if (rtype == -1) cycle
               call natashamcloane(isotot,kk,atot,aanf(i),ztot,zznf(i),nbranch(kk,rtype),amin,aareac)
               exit
            end do
         end do

         ! weak rates from netgen (mix of Takahashi & Yokoi 1987 and Goriely 1999)
         call netgen_interpolate_rates(t9, rho, ye, beta)

         do kk = 1, isotot
            if (atot(kk) > 226 .or. ztot(kk) > 85) cycle ! max A and Z for netgen weak rates

            do rtype = i_bm, i_ec ! beta-decay and electron capture
               i = netgen_beta_ntrans(atot(kk),ztot(kk),1,rtype)
               if (i <= 0) cycle
               vbranch(kk,rtype) = beta % rates(i)
               label(kk,rtype) = 'NETB1'
            end do
         end do

         ! weak targets from netgen:
         do i = 1, beta % ireac
            do kk = 1, isotot
               ! search for the isotope is in our network
               if (ztot(kk) /= beta%zt(i) .or. atot(kk) /= beta%at(i)) cycle

               ! if we find it, query netgen for the type of reaction this is
               rtype = netgen_rtype(beta, i)

               select case(rtype)
               case(i_bm,i_ec)
                  call natashamcloane(isotot,kk,atot,beta%af(i),ztot,beta%zf(i),nbranch(kk,rtype),amin,aareac)
               end select
               exit
            end do
         end do

         ! beta decay rates from fuller & fowler 1985 (light isotopes, up to Fe)
         ! do not use for decays!
         if (.not. decay) then
            call fuller_interpolate_rates(t9, rho, ye)

            ! assign rates at their position according to reaction type.
            do kk = 1, isotot
               if (atot(kk) > 60 .or. ztot(kk) > 30) cycle ! max a and z for ffn

               do rtype = i_bm, i_ec ! beta decay and electron capture
                  i = fuller_ntrans(atot(kk),ztot(kk),1,rtype)
                  if (i <= 0) cycle
                  vbranch(kk,rtype) = fuller_rates(i)
                  label(kk,rtype) = 'FFW85'
                  call arrow_check(index_raw,ztot,atot,eltot,kk,rtype,nbranch(kk,rtype))
               end do
            end do
         end if ! not decay

         ! weak rates from oda et al. 1994 (netgen).
         ! do not use for decays!
         if (.not. decay) then
            call netgen_interpolate_rates(t9, rho, ye, oda)

            do kk = 1, isotot
               if (atot(kk) > 39 .or. ztot(kk) > 20) cycle ! max a and z in oda 94 tables

               do rtype = i_bm, i_ec ! beta-decay and electron capture
                  i = netgen_oda_ntrans(atot(kk),ztot(kk),1,rtype)
                  if (i <= 0) cycle
                  vbranch(kk,rtype) = oda % rates(i)
                  label(kk,rtype) = 'ODA94'
                  call arrow_check(index_raw,ztot,atot,eltot,kk,rtype,nbranch(kk,rtype))
               end do
            end do
         end if ! not decay

         ! weak rates from Jones, Bertolli & Johnson 2016/2017.
         select case( jbj_mode )
         case(1)
            continue
         case(2,3,4)
            call jbj_interpolate_rates(ye, rho, t9)
            do kk = 1, isotot
               if (atot(kk) > jbjmaxA .or. ztot(kk) > jbjmaxZ) cycle
               if (atot(kk) < jbjminA .or. ztot(kk) < jbjminZ) cycle

               do rtype = i_bm, i_ec ! beta-decay and electron capture

                  i = jbj_ntrans(atot(kk),ztot(kk),1,rtype)
                  if (i <= 0) cycle
                  if (jbj_mode == 2 .and. label(kk,rtype) /= 'ODA94') cycle
                  vbranch(kk,rtype) = jbj_rate(i)
                  label(kk,rtype) = 'JBJ16'
                  call arrow_check(index_raw,ztot,atot,eltot,kk,rtype,nbranch(kk,rtype))
               end do
            end do
         case default
            stop "invalid jbj_mode"
         end select

         ! weak rates from netgen/langanke & martinez-pinedo 2000.
         ! do not use for decay calculations
         if (.not. decay) then
            call netgen_interpolate_rates(t9, rho, ye, lmp)

            do kk = 1, isotot
               if (atot(kk) > 65 .or. ztot(kk) > 32) cycle

               do rtype = i_bm, i_ec ! beta-decay and electron capture
                  i = netgen_lmp_ntrans(atot(kk),ztot(kk),1,rtype)
                  if (i <= 0) cycle
                  vbranch(kk,rtype) = lmp % rates(i)
                  label(kk,rtype) = 'LMP00'
                  call arrow_check(index_raw,ztot,atot,eltot,kk,rtype,nbranch(kk,rtype))
               end do
            end do
         end if ! not decay

         ! weak rates from Nabi, K-K & Volker 2004
         select case( nkk_mode )
         case(1)
            continue
         case(2,3)
            call nkk_interpolate_rates(rho, ye, t9)
            do kk = 1, isotot
               if (atot(kk) > nkkmaxA .or. ztot(kk) > nkkmaxZ) cycle
               if (atot(kk) < nkkminA .or. ztot(kk) < nkkminZ) cycle

               do rtype = i_bm, i_ec ! beta-decay and electron capture
                  i = nkk_ntrans(atot(kk),ztot(kk),1,rtype)
                  if (i <= 0) cycle
                  if (nkk_mode == 2 .and. label(kk,rtype) /= 'LMP00') then
                     ! only use NKK to replace LMP00 rates
                     cycle
                  end if
                  if (nkk_mode == 3) then
                     ! only use NKK where FFN, LMP and Oda94 are not available
                     if (label(kk,rtype) == "FFW85") cycle
                     if (label(kk,rtype) == "LMP00") cycle
                     if (label(kk,rtype) == "ODA94") cycle
                  end if
                  vbranch(kk,rtype) = nkk_rate(i)
                  label(kk,rtype) = 'NKK04'
                  call arrow_check(index_raw,ztot,atot,eltot,kk,rtype,nbranch(kk,rtype))
               end do
            end do
         case default
            stop "invalid choice for nkk_mode"
         end select

         ! nacre rates and targets from netgen.
         call netgen_interpolate_rates(t9, rho, ye, nacre)

         do i = 1, nacre % ireac
            do kk = 1,isotot
               ! search for the isotope is in our network
               if (ztot(kk) /= nacre%zt(i) .or. atot(kk) /= nacre%at(i)) cycle

               ! if we find it, query netgen for the type of reaction this is
               rtype = netgen_rtype(nacre, i)

               select case(rtype)
               case(i_pg)
                  aareac = 1
               case(i_pn)
                  aareac = 0
               case(i_pa)
                  aareac = 0
               case(i_ag)
                  aareac = 1
               case(i_an)
                  aareac = 0
               case(i_ap)
                  aareac = 0
               case default
                  cycle
               end select

               call natashamcloane(isotot,kk,atot,nacre%af(i),ztot,nacre%zf(i),nbranch(kk,rtype),amin,aareac)
               vbranch(kk,rtype) = nacre % rates(i)
               label  (kk,rtype) = nacre_name
               exit
            end do
         end do

         ! iliadis2001 proton captures and targets from netgen
         call netgen_interpolate_rates(t9, rho, ye, illi)

         do i = 1, illi % ireac
            do kk = 1, isotot
               ! search for the isotope is in our network
               if (ztot(kk) /= illi%zt(i) .or. atot(kk) /= illi%at(i)) cycle

               ! if we find it, query netgen for the type of reaction this is
               rtype = netgen_rtype(illi, i)

               select case(rtype)
               case(i_pg)
                  aareac = 1
               case(i_pa)
                  aareac = 0
               case default
                  cycle
               end select

               call natashamcloane(isotot,kk,atot,illi%af(i),ztot,illi%zf(i),nbranch(kk,rtype),amin,aareac)
               vbranch(kk,rtype) = illi % rates(i)
               label  (kk,rtype) = "ILI01"
               exit
            end do
         end do

         ! kadonis database
         call kadonis_interpolate_rates(t9,rho)

         ! define rate and target from kadonis
         do i = 1, isotot
            do j = 1, nkad
               if (zkad(j) /= ztot(i) .or. akad(j) /= atot(i)) cycle
               vbranch(i,i_ng) = vkad(j)
               label  (i,i_ng)   = 'KADON'
               exit
            end do
         end do

         ! other reactions module
         if (use_other_nuc) then
            call calculate_other_nuc(t9,rho,ye)

            do i = 1, isotot
               do j = 1, num_other
                  if (other_z(j) /= ztot(i) .or. other_a(j) /= atot(i)) cycle
                  vbranch(i,other_rtype(j)) = other_rate(j)
                  label  (i,other_rtype(j)) = other_label(j)
               end do
            end do
         end if

         ! write isotopes names into speneu array, that will become zis
         do kk = 1, isotot
            write( speneu(kk),'(A2,I3)') eltot(kk), atot(kk)
            if ( eltot(kk) == 'nn') speneu(kk) = 'NEUT '
            if ( eltot(kk) == 'H '.and. atot(kk) == 1) speneu(kk) = 'PROT '
         end do

   end subroutine networkI



   !> this subroutine's job is to get the reaction rates for each
   !> reaction type for each species (vbranch) and return it. vbranch is
   !> of dimension (nre), so we get one rate per reaction. Sounds
   !> sensible, except for in the networkI routine there is also an
   !> array vbranch with dimension (nsp,<number of reaction types>) that
   !> holds the same information but in a different format, which causes
   !> a bit of confusion
   subroutine networkII ( ye, t9, rho, vbranch, NRNC, ilabb )
         use array_sizes
         use rates, only: v, locate_reaclib_in_vbranch, &
               locate_beta_netgen_in_vbranch, &
               locate_fuller_in_vbranch, locate_oda_netgen_in_vbranch, &
               locate_lmp_netgen_in_vbranch, locate_nacre_netgen_in_vbranch, &
               locate_illi_netgen_in_vbranch, locate_kadonis_in_vbranch, &
               locate_jbj_in_vbranch, locate_nkk_in_vbranch, &
               locate_other_in_vbranch
         use reaction_info, only: lab, labb
         use netgen, only: netgen_interpolate_rates, t9_netgen_limit, &
               beta, illi, oda, lmp, nacre
         use fuller, only: fuller_interpolate_rates, fuller_rates
         use physics_knobs
         use frame_knobs
         use utils
         use reaclib, only: aan, zzn, aanf, zznf, reaclib_ntrans, nireac, jtotdim, reaclib_name, &
               reaclib_wrapper, rv
         use kadonis, only: nkad, vkad, kadonis_interpolate_rates
         use jbj16, only: jbj_rate, jbj_interpolate_rates
         use nkk04, only: nkk_rate, nkk_interpolate_rates
         use other_nuc, only: calculate_other_nuc, other_rate
         use nuc_data
         use constants
         implicit none

         integer, dimension(nre) :: ilabb(nre)
         integer  :: NRNC, i, j, kk, ireac

         real(r8) :: t9, ye, rho, tti1, tti2, mrg_time
         real(r8) :: vbranch(nre)

         ! profiling
         real(r8) :: rlb_time, bta_time, ffn_time, lmp_time, oda_time, &
               kds_time, jbj_time, nkk_time, nac_time, ill_time, other_nuc_time

         ! reaclib

         tti1 = wallclocktime()
         call reaclib_wrapper(t9, rho, ye)
         rlb_time = wallclocktime() - tti1

         ! beta decay rates from netgen

         tti1 = wallclocktime()
         call netgen_interpolate_rates(t9, rho, ye, beta)
         bta_time = wallclocktime() - tti1

         ! beta decay rates from fuller & fowler 1985 (light isotopes, up to Fe)

         tti1 = wallclocktime()
         call fuller_interpolate_rates(t9, rho, ye)
         ffn_time = wallclocktime() - tti1

         ! oda et al. 1994

         tti1 = wallclocktime()
         call netgen_interpolate_rates(t9, rho, ye, oda)
         oda_time = wallclocktime() - tti1

         ! Jones, Bertolli & Johnson 2016/2017.

         if (jbj_mode /= 1) then
            tti1 = wallclocktime()
            call jbj_interpolate_rates(ye, rho, t9)
            jbj_time = wallclocktime() - tti1
         end if

         ! Langanke & Martinez-Pinedo 2000

         tti1 = wallclocktime()
         call netgen_interpolate_rates(t9, rho, ye, lmp)
         lmp_time = wallclocktime() - tti1

         ! Nabi, K-K & Volker 2004

         if (nkk_mode /= 1) then
            tti1 = wallclocktime()
            call nkk_interpolate_rates(rho, ye, t9)
            nkk_time = wallclocktime() - tti1
         end if

         ! other_nuc module

         if (use_other_nuc) then
            tti1 = wallclocktime()
            call calculate_other_nuc(t9, rho, ye)
            other_nuc_time = wallclocktime() - tti1
         end if

         ! below t9_netgen_limit netgen tables are out of range.
         ! then, no proton captures or neutron captures are activated

         if (t9 > t9_netgen_limit) then

            ! nacre rates from netgen

            tti1 = wallclocktime()
            call netgen_interpolate_rates(t9, rho, ye, nacre)
            nac_time = wallclocktime() - tti1

            ! iliadis2001 from netgen;  (p,g ),(p,a) reactions included

            tti1 = wallclocktime()
            call netgen_interpolate_rates(t9, rho, ye, illi)
            tti2 = wallclocktime()
            ill_time = wallclocktime() - tti1

            ! kadonis

            tti1 = wallclocktime()
            call kadonis_interpolate_rates(t9,rho)
            tti2 = wallclocktime()
            kds_time = wallclocktime() - tti1
         else
            nacre % rates(:) = ZERO
            illi % rates (:) = ZERO
            vkad         (:) = ZERO
         end if ! t9 > t9_netgen_limit


         ! during initialisation, the positions (*ntrans(ant(kk),znt(kk)+1,1,ilabb(kk))) are saved, to
         ! prevent 4d array access during the main cycle loop.
         tti1 = wallclocktime()

         do kk = 1, nrnc
            select case(lab(kk))
            case('BASEL', 'JINAR', 'JINAC', 'JINAV')
               i = locate_reaclib_in_vbranch(kk)
               vbranch(kk) = rv(i)
            case('NETB1')
               i = locate_beta_netgen_in_vbranch(kk)
               vbranch(kk) = beta % rates(i)
            case('FFW85')
               i = locate_fuller_in_vbranch(kk)
               vbranch(kk) = fuller_rates(i)
            case('ODA94')
               i = locate_oda_netgen_in_vbranch(kk)
               vbranch(kk) = oda % rates(i)
            case('LMP00')
               i = locate_lmp_netgen_in_vbranch(kk)
               vbranch(kk) = lmp % rates(i)
            case('NACRL','NACRR','NACRU')
               i = locate_nacre_netgen_in_vbranch(kk)
               vbranch(kk) = nacre % rates(i)
            case('ILI01')
               i = locate_illi_netgen_in_vbranch(kk)
               vbranch(kk) = illi % rates(i)
            case('KADON')
               i = locate_kadonis_in_vbranch(kk)
               vbranch(kk) = vkad(i)
            case('JBJ16')
               i = locate_jbj_in_vbranch(kk)
               vbranch(kk) = jbj_rate(i)
            case('NKK04')
               i = locate_nkk_in_vbranch(kk)
               vbranch(kk) = nkk_rate(i)
            case('OOOOO','VITAL','MPG06','NTRNO','RVRSE')
               ! blank, vital, a-decay, nautrino or reverse reaction
               continue
            case default
               i = locate_other_in_vbranch(kk)
               vbranch(kk) = other_rate(i)
            end select
         end do

         mrg_time = wallclocktime() - tti1

         ! now we have all the vbranch (reacion rates) at the right temperature and density

#ifdef PROFILE
         write (*,fmtprofile) &
               'networkII/rlb cpu_time/s = ',   rlb_time
         write (*,fmtprofile) &
               'networkII/bta cpu_time/s = ',   bta_time
         write (*,fmtprofile) &
               'networkII/ffn cpu_time/s = ',   ffn_time
         write (*,fmtprofile) &
               'networkII/lmp cpu_time/s = ',   lmp_time
         write (*,fmtprofile) &
               'networkII/oda cpu_time/s = ',   oda_time
         write (*,fmtprofile) &
               'networkII/nac cpu_time/s = ',   nac_time
         write (*,fmtprofile) &
               'networkII/ill cpu_time/s = ',   ill_time
         write (*,fmtprofile) &
               'networkII/kds cpu_time/s = ',   kds_time
         write (*,fmtprofile) &
               'networkII/jbj cpu_time/s = ',   jbj_time
         write (*,fmtprofile) &
               'networkII/nkk cpu_time/s = ',   nkk_time
         write (*,fmtprofile) &
               'networkII/other cpu_time/s = ', other_nuc_time
         write (*,fmtprofile) &
               'networkII mrg cpu_time/s = ',   mrg_time
#endif
   end subroutine networkII


   subroutine rnetw2007( ye, qi, an, zn, nvar, nvrel, rho, t9, yps )
         ! The whole network setup is controlled in this routine. It is
         ! only ever called during setup (and only then if a network has
         ! not been defined already in networksetup.txt by the user).
         ! So, the variable istart should always be 0, meaning "init"
         use array_sizes
         use nuc_data
         use rates
         use physics_knobs
         use frame_knobs
         use utils
         use reaction_info, only: num_rtypes, i_ng, i_gn, i_np, i_na, i_pg, i_gp, i_pn, i_pa, &
              ! RJS 23/11/19 -- added i_bp
               i_ag, i_ga, i_an, i_ap, i_bm, i_ec, i_bn, i_bp, i_ba, i_ve, i_vp, i_vsp, i_vsn, i_vsa, rlabels
         use vital, only: vital_index, read_physics_input_data, vital_rates_derivs
         use networksetup, only: nvnc1, nrnc1, inumer, inumerr
         use alpha_decays
         use neutrinos
         use constants
         use communication
         use reaction_info, only: lab, labb, ilabb, bind_energy_diff, rfac
         implicit none

         ! nvcp Number isotopes in charged particle network
         ! nvnc ................ in neutron capture network
         ! NR.......... reactions ....
         integer :: k,L,i,kk,j, & !< various counters
               ispecies, & !< species index
               rtype, & !< reaction type, e.g. (p,g), (a,p), etc
               nvrel,nvar,ngir,ngr,ngis,ngs,nvnc,nrnc
         integer, allocatable :: nbranch(:,:)
         real(r8) :: an(nsp), zn(nsp), qi(nre), ye, rho, t9, yps(nsp)
         real(r8), allocatable :: vbranch(:,:)
         character(len=5) speneu(nsp)
         character(len=5), allocatable :: label(:,:)

         ! only considerisotope, considerreaction and V from new network are
         ! reset to zero.
         real(r8) :: tti1, tti2
         integer :: &
               iprot, &  !< proton index
               ineut, &  !< neutron index
               ihe4,  &  !< he4 index
               igamma    !< gamma index

         ! initialisation
         allocate(nbranch(nsp,num_rtypes), vbranch(nsp,num_rtypes), label(nsp,num_rtypes))
         considerisotope  ( nvcp+1: ) = .false.
         considerreaction ( nrcp+1: ) = .false.
         nvnc = 0; nrnc = 0; v(:) = ZERO

         ! read species and reaction data for VITAL from ppn_physics.input
         call read_physics_input_data(an, zn, qi, nrcp, nvcp)

         ! give names to the vital rates index, which can be used
         ! to substitute the vital rates with reaclib rates in networkII.
         call vital_index(nrcp)

         ! activate terrestrial decay for species in vital with half life lower than tbetamin.
         do i = 1, nvcp
            if (.not. considerisotope(i)) cycle
            if (t_half(i) < tbetamin) then
               do l = 1, nrcp
                  if (ispe(zis(i)) == k1(l)) then
                     if (ilabb(l) == i_bm .or. ilabb(l) == i_ec) then
                        considerreaction(l)= .true.
                        write(*,*) &
                              'terrestrial decay rate activated in vital'
                        write(*,*)'isotope:',zis(i)
                        exit
                     end if
                  end if
               end do
            end if
         end do

         qi(:) = ZERO
         iprot = ispe('PROT '); ineut = ispe('NEUT '); ihe4 = ispe('HE  4'); igamma = ispe('OOOOO')


         ! now fill isotope and reaction list with those from neutron network
         ! we do not replace already existing reactions in previous NW!
         ! but we may add some reactions that are already existant from use
         ! of nnn<0

         ! set counters to number of isotopes/reactions in "charged particle" network
         k     = nvcp ! new isotopes
         L     = nrcp ! new reactions
         nvnc  = nvcp

         if (nnn > 0) then
            !> nnn is the number of isotopes that should be included from the file isotopedatabase.txt.
            !> therefore, nnn > 0 indicates that we should call the subroutine networkI, which defines
            !> all reactions and isotopes that are not included in ppn_physics.input
            nbranch(:,:) = 0
            vbranch(:,:) = ZERO

            call networkI (ye, t9, rho, nbranch, vbranch, speneu, label)

            do i = 1, nnn
               ispecies = ispe(speneu(i))
               if (ispecies > nvcp .or. ispecies == -1) then
                  k       = k + 1
                  kk      = k
                  nvnc    = k

                  zis(kk) = speneu(i)

                  read(zis(kk)(3:5),'(F3.0)') AN(kk)

                  zn(kk)  = ztot(i)
                  considerisotope(kk) = .true.
               else
                  kk = ispecies
               endif
               do rtype = 1, num_rtypes
                  if (nbranch(i,rtype) == 0) cycle

                  L = L + 1
                  if (L > nre) stop 'ERROR in rnetw2007: L > nre'
                  considerreaction(L) = .true.

                  k1(L)    = kk
                  k7(L)    = nbranch(i,rtype)
                  v(L)     = vbranch(i,rtype)
                  lab(L)   = label  (i,rtype)
                  ilabb(L) = rtype
                  labb(L)  = rlabels(rtype)
                  k2(L)    = 1
                  k8(L)    = 1

                  select case(rtype)
                  case(i_ng)
                     k3(L) = ineut ; k4(L) = 1; k5(L) = igamma; k6(L) = 0
                  case(i_gn)
                     k3(L) = igamma; k4(L) = 0; k5(L) = ineut ; k6(L) = 1
                  case(i_np)
                     k3(L) = ineut ; k4(L) = 1; k5(L) = iprot ; k6(L) = 1
                  case(i_na)
                     k3(L) = ineut ; k4(L) = 1; k5(L) = ihe4  ; k6(L) = 1
                  case(i_pg)
                     k3(L) = iprot ; k4(L) = 1; k5(L) = igamma; k6(L) = 0
                  case(i_gp)
                     k3(L) = igamma; k4(L) = 0; k5(L) = iprot ; k6(L) = 1
                  case(i_pn)
                     k3(L) = iprot ; k4(L) = 1; k5(L) = ineut ; k6(L) = 1
                  case(i_pa)
                     k3(L) = iprot ; k4(L) = 1; k5(L) = ihe4  ; k6(L) = 1
                  case(i_ag)
                     k3(L) = ihe4  ; k4(L) = 1; k5(L) = igamma; k6(L) = 0
                  case(i_ga)
                     k3(L) = igamma; k4(L) = 0; k5(L) = ihe4  ; k6(L) = 1
                  case(i_an)
                     k3(L) = ihe4  ; k4(L) = 1; k5(L) = ineut ; k6(L) = 1
                  case(i_ap)
                     k3(L) = ihe4  ; k4(L) = 1; k5(L) = iprot ; k6(L) = 1
                  case(i_bm)
                     k3(L) = igamma; k4(L) = 0; k5(L) = igamma; k6(L) = 0
                  case(i_ec)
                     k3(L) = igamma; k4(L) = 0; k5(L) = igamma; k6(L) = 0
                  case(i_bn)
                     k3(L) = igamma; k4(L) = 0; k5(L) = ineut ; k6(L) = 1
                  ! RJS 23/11/19 -- added case for beta_delayed proton emission
                  case(i_bp)
                     k3(L) = igamma; k4(L) = 0; k5(L) = iprot ; k6(L) = 1
                  case default
                     stop 'unknown reaction type in rnetw2007'
                  end select
               end do
            enddo
         end if ! (nnn>0)

         nvnc = nvnc - nvcp
         NRNC = L - nrcp

         ! identify real k7:
         do l = nrcp + 1, nrcp + NRNC
            k7(l) = ispe(speneu(k7(l)))
            ! in reaclib sources reactions like 2p+00000 are written like
            ! 1p+1p, and their format is not correct for solver stuff.
            ! here below we correct this things.
            if (k1(l) == k3(l)) then
               k2(l) = k2(l) + k4(l)
               k4(l) = 0
               k3(l) = igamma
            end if
            ! the same is done for k5 and k7, after real k7 is defined, e.g. Be7(n,a)a
            if (k5(l) == k7(l)) then
               k8(l) = k8(l) + k6(l)
               k6(l) = 0
               k5(l) = igamma
            end if
         end do

         ! alpha decays are not included in reaclib, so they are out of
         ! the main network we have built so far...
         ! alpha decay by Magill, J., Pfennig, G. and Galy, J. 2006 added
         ! (this is Karlsruhe isotopic chart 2006 edition, revised 2007)
         iadcy = nrcp + nrnc

         ! when the subroutine natashamcloane was designed, I was not considering that
         ! the most favourable decay can be alpha-decay.... therefore, in case is missing
         ! an isotope, the network would assume beta decay and not alpha.
         ! Until this is not fixed, I just exclude all decays with no isotopes F involved.
         ! Also, before using considerisotope we do not check yet if p_iso and d_iso
         ! are actually in isotopedatabase.... here below we exclude reactions involving
         ! isotopes not included in isotopedatabase.txt or included as false.
         if (.not. adcy_masks_exist) call alpha_decays_create_masks()
         do i = 1, num_packed_adcy
            k1    (iadcy+i) = ispe(p_iso_packed(i))
            k2    (iadcy+i) = 1
            k3    (iadcy+i) = igamma
            k4    (iadcy+i) = 0
            k5    (iadcy+i) = ihe4
            k6    (iadcy+i) = 1
            k7    (iadcy+i) = ispe(d_iso_packed(i))
            k8    (iadcy+i) = 1
            v     (iadcy+i) = adcy_rate_packed(i)
            lab   (iadcy+i) = adcy_name
            labb  (iadcy+i) = '(b,a)'
            ilabb (iadcy+i) = i_ba ! alpha-decay
            rfac  (iadcy+i) = 1._r8

            bind_energy_diff(iadcy+i) = ZERO
            considerreaction(iadcy+i) = .true.
         end do

         ! updating reaction count for the alpha decays
         nrnc = nrnc + num_packed_adcy

         ! neutrino reactions
         if (do_neutrinos) then
            if (.not. nu_masks_exist) call nu_create_masks()
            i_nu = nrcp + nrnc
            k = 0
            do i = 1, num_nu_reac
               if (.not. nu_mask(i)) cycle
               k = k + 1
               k1    (i_nu+k) = ispe(epsi(nu_a(i),nu_z(i)))
               k2    (i_nu+k) = 1
               k3    (i_nu+k) = igamma
               k4    (i_nu+k) = 0
               select case(nu_type(i))
               case(i_nue)
                  ! electron neutrino capture
                  k5    (i_nu+k) = igamma
                  k6    (i_nu+k) = 0
                  k7    (i_nu+k) = ispe(epsi(nu_a(i),nu_z(i)+1))
                  k8    (i_nu+k) = 1
                  labb  (i_nu+k) = '(v,-)'
                  ilabb (i_nu+k) = i_ve
               case(i_nup)
                  ! antielectron neutrino capture
                  k5    (i_nu+k) = igamma
                  k6    (i_nu+k) = 0
                  k7    (i_nu+k) = ispe(epsi(nu_a(i),nu_z(i)-1))
                  k8    (i_nu+k) = 1
                  labb  (i_nu+k) = '(v,+)'
                  ilabb (i_nu+k) = i_vp
               case(i_nusp)
                  ! proton spallation
                  k5    (i_nu+k) = iprot
                  k6    (i_nu+k) = 1
                  k7    (i_nu+k) = ispe(epsi(nu_a(i)-1,nu_z(i)-1))
                  k8    (i_nu+k) = 1
                  labb  (i_nu+k) = '(v,p)'
                  ilabb (i_nu+k) = i_vsp
               case(i_nusn)
                  ! neutron spallation
                  k5    (i_nu+k) = iprot
                  k6    (i_nu+k) = 1
                  k7    (i_nu+k) = ispe(epsi(nu_a(i)-1,nu_z(i)))
                  k8    (i_nu+k) = 1
                  labb  (i_nu+k) = '(v,n)'
                  ilabb (i_nu+k) = i_vsn
               case(i_nusa)
                  ! alpha spallation
                  k5    (i_nu+k) = ihe4
                  k6    (i_nu+k) = 1
                  k7    (i_nu+k) = ispe(epsi(nu_a(i)-2,nu_z(i)-2))
                  k8    (i_nu+k) = 1
                  labb  (i_nu+k) = '(v,a)'
                  ilabb (i_nu+k) = i_vsa
               end select
               lab   (i_nu+k) = nu_name
               v     (i_nu+k) = v_nu(i)
               rfac  (i_nu+k) = 1._r8

               bind_energy_diff(i_nu+k) = ZERO
               considerreaction(i_nu+k) = .true.
            end do
            ! updating reaction count for the neutrino reactions
            nrnc = nrnc + num_packed_nu
         end if

         ! exclude all reactions from rnetw2007 network in common with vital subroutine
         do i = 1, nrcp
            do j = nrcp + 1, nrcp + nrnc
               if (zis(k1(i)) == zis(k1(j)) .and. k2(i) == k2(j) .and.  k3(i) == k3(j) .and. k4(i) == k4(j)) then
                  if(k8(i)==k8(j).and.k6(i) == k6(j).and.zis(k5(i))==zis(k5(j)) .and. considerreaction(i)) then
                     considerreaction(j) = .false.
                  end if
               end if
            end do
         end do

         ! exclude all reactions from rnetw2007 network that include species
         ! that are false in vital network (e.g., NEUT).
         do i = 1, nvcp
            if (.not. considerisotope(i) .and. ispe(zis(i)) /= igamma) then
               if (master) write(*,*)'vital: false species = ', zis(i)
               do L = 1, nrcp + NRNC
                  if (ispe(zis(i)) == k1(L) .or. ispe(zis(i)) == k3(L)) then
                     considerreaction(L) = .false.
                  end if
                  if (ispe(zis(i)) == k5(L) .or. ispe(zis(i)) == k7(L)) then
                     considerreaction(L) = .false.
                  end if
               end do
            end if
         end do

         ! insert here part to build isomer frame:
         ! 1) I build isomer reactions as twin of ground state reactions.
         ! 2) proper beta decay rate is given to ground and isomer.
         ! ngs = number of ground species;
         ! ngis = number of ground + isomer species;
         ! ngr = number of ground species reactions;
         ! ngir = number of ground + isomer species reactions.

         ngs = nvcp + nvnc
         ngr = nrcp + nrnc

         call isomerI( t9, an, zn, lab, labb, ilabb, ngs, ngr, inumer)
         !      call isomerII( t9, an, zn, lab, labb, ilabb, ngs, ngr, inumer)

         ngis     = ngs ! total number of species including isomers?
         ngir     = ngr ! total number of reactions including isomers?
         INUMERR  = ngir - nrcp - NRNC ! number of isomer reactions?

         ! compute NVAR and NVREL including isomers:
         nvar     = count(considerisotope)  ! number of active (.true.) species in network
         nvrel    = count(considerreaction) ! number of active (.true.) reactions in network

         if (master) then
            if (iolevel >= 3) print *, 'rnetw2007: nvar =', nvar,' and nvrel =', nvrel
         end if

         ! creation of ant, znt. used in rnetw2008 to merge rates in ppn network.
         do i = 1, ngis
            do l = 1, ngir
               if (zis(i) == zis(k1(l))) then
                  ant(l) = int(an(i))
                  znt(l) = int(zn(i))
               end if
               if (zis(i) == zis(k3(l))) then
                  anf(l) = int(an(i))
                  znf(l) = int(zn(i))
               end if
            end do
         end do

         ! creation of index array niso
         ! third argument is 1 (ground) or 2 (isomeric).
         ! It will account of the ground/isomeric state.
         niso(:,:,:) = 0
         do i = 1, nvcp + nvnc
            isomeric_state(i) = 1
            if (considerisotope(i)) then
               niso(nint(an(i)),nint(zn(i)),isomeric_state(i)) = i
            end if
         end do
         if ( inumer > 0 ) then
            do i = nvcp + nvnc + 1, ngis
               isomeric_state(i) = 2
               if (considerisotope(i))then
                  niso(nint(an(i)),nint(zn(i)),isomeric_state(i)) = i
               end if
            end do
         end if

         ! hash the locations of each rate source in the main rate array
         call rates_hash_locations_for_merge(ngir, ilabb)

         ! calculate the binding energy variation associated to each reaction rate.
         call get_binding_energy_diff(an,zn,inumer,ngis)

         ! write out what has been selected
         nvnc1 = ngis ; nrnc1 = ngir
         ! set rate factors for "non-charged-particle" (just non-VITAL I think)
         rfac(nrcp+1:ngir) = ONE

         ! get vital rates so that the entries in networksetup.txt are initialised
         call vital_rates_derivs(rho_nw_ini, t9_nw_ini, yps)

         deallocate(nbranch); deallocate(vbranch); deallocate(label)

   end subroutine rnetw2007



   subroutine rnetw2008( ye, an, zn, nvar, nvrel, rho, t9, x, nu )
         use communication
         use array_sizes
         use nuc_data
         use rates
         use physics_knobs
         use frame_knobs
         use utils
         use vital, only: vital_index, vital_rates_derivs
         use networksetup, only: nvnc1, nrnc1, inumer, inumerr, &
               write_networksetup2
         use alpha_decays
         use neutrinos
         use constants
         use reaction_info, only: lab, labb, ilabb, bind_energy_diff, rfac
         implicit none

         integer i, j, k, l, nvrel, nvar
         real(r8) :: an(nsp), zn(nsp), ye, tti1, tti2
         integer :: &
               nvnc, & !< number isotopes in "neutron capture network"
               nrnc    !< number reactions in "neutron capture network"

         ! (nvcp: Number isotopes in charged particle network)
         ! (nrcp: Number reactions in charged particle network)
         real(r8) :: vbranch(nre), rho, t9, x(nsp), vit_time
         type(neutrino) :: nu
            
         integer nrcpold
         common/old/nrcpold

         ! initialise
         v (:) = ZERO

         ! get vital rates

         tti1 = wallclocktime()
         call vital_rates_derivs(rho, t9, x)
         vit_time = wallclocktime() - tti1
#ifdef profile
         write (*,fmtprofile) 'vital cpu_time/s = ',   vit_time
#endif

         nvnc = nvnc1
         nrnc = nrnc1

         if (nnn > 0) then
            ! get rates beyond VITAL
            call networkII (ye, t9, rho, vbranch, nrnc-inumerr, ilabb)
            v(nrcpold+1:nre) = vbranch(nrcpold+1:nre)
         end if

         ! get alpha-decay rates and insert them

         iadcy = nrnc - inumerr - num_packed_adcy - num_packed_nu
         if (.not. adcy_masks_exist) call alpha_decays_create_masks()
         v(iadcy+1:iadcy+num_packed_adcy) = adcy_rate_packed(1:num_packed_adcy)

         ! get neutrino rates and appand them
         if (do_neutrinos) then
            i_nu = nrnc - inumerr - num_packed_nu
            if (.not. nu_masks_exist) call nu_create_masks()
            call calculate_neutrino_rates(t9, rho, ye, nu)
            v(i_nu+1:i_nu+num_packed_nu) = v_nu_packed(:)
         end if

         ! insert here part to build isomer frame:
         !    1) I build isomer reactions as twin of ground state reactions,
         !       if not thermalized.
         !    2) proper beta decay rate is given to ground and isomer.

         nvar = count(considerisotope)
         call isomerI( t9, an, zn, lab, labb, ilabb, nvnc-inumer, nrnc-inumerr, inumer )
         !      call isomerII( t9, an, zn, lab, labb, ilabb, nvnc-inumer, nrnc-inumerr, inumer )

         ! recompute new NVAR and new NVREL:

         nvar = count(considerisotope)
         nvrel = count(considerreaction)

         ! reporting for high IO levels:

         if (iolevel >= 4) print *, "Physics: nvar after isomerI", nvar
         if (iolevel >= 3) print *, 'rnetw2008: nvar =',nvar,' and nvrel ='

         if (ininet == 3 .and. istart == 0 .and. master) then
            print *, "TODO: The ininet=3 option is currently not working."
            print *, "Reaction rate factors can be implemented via the physics namelist ", &
                  "see physics_knobs.F90."
            !stop
            !call write_networksetup2(t9, rho, nvar, nrcp, nvrel, an, zn)
            !call write_networksetup2(t9, rho, nvnc, nrnc, nvrel, an, zn)
         end if

   end subroutine rnetw2008



   subroutine isomerI( t9, an, zn, lab, labb, ilabb, nums, numr, n_isomers_in_network )
         use utils, only: r8
         use array_sizes, only: nsp, nre
         use frame_knobs, only : istart
         use isomers, only: isomers_init, isomers_calculate_rates, &
               num_active_isomers

         implicit none

         real(r8) :: t9, an(nsp), zn(nsp)
         integer :: nums, numr, ilabb(nre), n_isomers_in_network
         character(len=5) :: lab(nre), labb(nre)

         if (istart == 0) then
            ! TODO: these lab, ilab and ilabb need more transparent name
            ! and probably don't need to be passed here (should be taken
            ! isomer module directory from wherever they live)
            call isomers_init(an, zn, nums, numr, lab, labb, ilabb)

            n_isomers_in_network = num_active_isomers
         end if

         if (n_isomers_in_network == 0) return

         call isomers_calculate_rates(t9, numr, lab, labb, ilabb)
   end subroutine isomerI


   !> translate reaction type (bm = 13, ec = 14) into (bm = 1, ec = -1)
   integer(i4) function direction(rtype)
         integer(i4) :: rtype
         direction = -(2 * rtype - 27)
   end function direction

   !> this subroutine gives the targets (product nucleus) for weak reactions. If the target 
   !> does not exist, it will assume that the target would have instantaneously beta-decayed
   !> and set the target to the second beta decay product, e.g. 89Sr --> 89Zr. If the second
   !> beta-decay product does not exist, it will go to third, and so on, for 4 links in the 
   !> beta-decay chain. If there is still no valid target, the reaction will be excluded
   !> from the network and a warning to that effect will be written to the terminal
   subroutine arrow_check(index_raw,ztot,atot,eltot,ind,rtype,produced_species)
         use array_sizes
         use frame_knobs
         use physics_knobs, only: index_reaclib
         implicit none

         integer :: atot(nsp), ztot(nsp), &
               rtype, & !< reaction type; see reaction_info module for list of types
               produced_species,ind,idirection_of_decay

#if pIDX_RCLB == 3
         integer :: index_raw(0:iCfdim,i325dim) !> index for weak reactions full
#else
         integer :: index_raw(0:iAtdim,i282dim) !> index for weak reactions full
#endif
         character(len=2) :: eltot(nsp)

         idirection_of_decay = direction(rtype)

         if (idirection_of_decay == 1) then
            if (produced_species == 0) produced_species = index_raw(ztot(ind)+1,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)+2,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)+3,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)+4,atot(ind))
            if (produced_species == 0) then
               if (iolevel .ge. 2) then
                  write(*,*)'warning: target of ',eltot(ind),atot(ind), &
                     'b- is not included. I am skipping the decay. check the network'
               end if
            end if
         else if(idirection_of_decay == -1)then
            if (produced_species == 0) produced_species = index_raw(ztot(ind)-1,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)-2,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)-3,atot(ind))
            if (produced_species == 0) produced_species = index_raw(ztot(ind)-4,atot(ind))
            if (produced_species == 0) then
               if (iolevel >= 2) then
                  write(*,*)'warning: target of ',eltot(ind),atot(ind), &
                     'b+ is not included. I am skipping the decay. check the network'
               end if
            end if
         end if
      end subroutine arrow_check



   subroutine get_binding_energy_diff(an,zn,inumer,ngis)
         use array_sizes
         use nuc_data
         use rates
         use reaclib, only: anumpart, znumpart, istate, spinpart, partnum, exmass
         use constants, only: ev2erg, mev2erg, avogadro, deltap, deltan, ZERO
         use utils, only: r8
         use reaction_info, only: bind_energy_diff
         implicit none

         integer i,j

         real(r8) :: an(nsp),zn(nsp),bind_energy(0:nsp),excess_mass(nsp)
         ! *** to get exmass data
         ! *** for also isomers
         integer ion_dum, inumer, ngis

         excess_mass = ZERO ; bind_energy = ZERO ; bind_energy_diff = ZERO
         j = 0 ; ion_dum = 0
         do i=1,nnpartdim
            if (istate(i).ge.1)then
               if(istate(i) == 1)ion_dum=1
               if(istate(i) == 2)ion_dum=1
               if(istate(i) == 3)ion_dum=1 !2
               ! I assume that thermlaized/gr/isomeric states have the same binding energies.
               ! in case there are isomers in these conditions. Marco Sep 2010
               j = niso(anumpart(i),znumpart(i),ion_dum)
               if(j.gt.0)then
                  excess_mass(j) = exmass(i)
               end if
            end if
         end do
         ! for now I do not have the excess mass for isomers.
         if (inumer .gt. 0)then
            do i=ngis-inumer+1,ngis
               excess_mass(i) = excess_mass(niso(nint(an(i)),nint(zn(i)),1))
            end do
         end if ! now I have them. :)

         do i=1,nsp
            if(considerisotope(i))then
               ! this provides really similar with what Frank has
               bind_energy(i) = zn(i)*deltap + (an(i) - zn(i))*deltan - excess_mass(i)
            end if
         end do
         ! I set binding energy for neut and prot = 0.
         bind_energy(ispe('NEUT ')) = ZERO
         bind_energy(ispe('PROT ')) = ZERO

         ! *** here I assign for a given reaction the binding energy variation,
         ! *** given the species that are involved.

         do i=1,nre
            if(considerreaction(i))then
               if (k7(i) == -1) then
                  print *, i, k1(i), k3(i), k5(i), k7(i)
                  stop
               end if
               bind_energy_diff(i) = &
                     -( &
                     k2(i)*bind_energy(niso(nint(an(k1(i))),nint(zn(k1(i))),isomeric_state(k1(i)))) + &
                     k4(i)*bind_energy(niso(nint(an(k3(i))),nint(zn(k3(i))),isomeric_state(k3(i)))) &
                     ) + &
                     k6(i)*bind_energy(niso(nint(an(k5(i))),nint(zn(k5(i))),isomeric_state(k5(i)))) + &
                     k8(i)*bind_energy(niso(nint(an(k7(i))),nint(zn(k7(i))),isomeric_state(k7(i))))
            end if
         end do

         bind_energy_diff = bind_energy_diff * mev2erg * avogadro
   end subroutine get_binding_energy_diff


end module physics
