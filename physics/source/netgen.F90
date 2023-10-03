module netgen
   use array_sizes, only : nre1, ngrid
   use utils, only: r8, netgen_data, irate, idrdt, idrdd, bsearch_r8
   use constants
   use communication

   implicit none
   private

   type(netgen_data), target, public :: &
         beta, oda, lmp, nacre, illi

   ! ^_^ for merging into ppn network
   integer, public :: &
         netgen_beta_ntrans (226, 0:85, 1, 13:14), &
         netgen_nacre_ntrans( 31, 0:15, 1,  4:13), &
         netgen_illi_ntrans ( 41, 0:21, 1,  5:8 ), &
         netgen_oda_ntrans  ( 39, 0:20, 1, 13:14), &
         netgen_lmp_ntrans  ( 65, 0:32, 1, 13:14)

   ! ^_^ file handles and file names for netgen rate files
   integer, parameter :: &
         netgen_beta_fh = 931, netgen_oda_fh = 932, netgen_lmp_fh = 933, netgen_nacre_fh = 934, &
         netgen_illi_fh = 935

   character(len=256), parameter :: &
         netgen_beta_fn     = '../NPDATA/netgen/netgen_betarecom_log_30_d_2012.txt', &
         netgen_oda_fn      = '../NPDATA/netgen/netgen_oda94_log_30_d2.txt', &
         netgen_lmp_fn      = '../NPDATA/netgen/netgen_lmp00_log_30_d.txt', &
         netgen_illi_fn     = '../NPDATA/netgen/netgen_iliadis2001_log_100.txt', &
         netgen_nacre1_fn   = '../NPDATA/netgen/netgen_nacre_log_100_low.txt', &
         netgen_nacre2_fn   = '../NPDATA/netgen/netgen_nacre_submit_rec.txt', &
         netgen_nacre3_fn   = '../NPDATA/netgen/netgen_nacre_log_100_up.txt'

   character(len=256)   :: nacrefn
   character(len=5)     :: nacre_name
   integer, parameter   :: iwhatnacre = 2

   ! tables in netgen are not given below 10^6 K, so set rate to 0 below
   ! shouldn't they be clipped to the rate at 1e6 K instead of set to 0?? Anyway, we don't to
   ! nucleosynthesis below 1e7 K
   real(r8), parameter, public  :: t9_netgen_limit = 0.001_r8

   public :: ngrid, nacre_name, iwhatnacre, netgen_init, netgen_data, netgen_interpolate_rates, &
         netgen_rtype


contains



   subroutine netgen_init()

         ! ^_^ open netgen files and read them into the netgen data structure for
         !     each file/source of data

         netgen_beta_ntrans(:,:,:,:)   = 0
         netgen_nacre_ntrans(:,:,:,:)  = 0
         netgen_illi_ntrans(:,:,:,:)   = 0
         netgen_oda_ntrans(:,:,:,:)    = 0
         netgen_lmp_ntrans(:,:,:,:)    = 0

         ! set names
         select case(iwhatnacre)
         case(1)
            nacrefn = netgen_nacre1_fn
            nacre_name = 'NACRL'
         case(2)
            nacrefn = netgen_nacre2_fn
            nacre_name = 'NACRR'
         case(3)
            nacrefn = netgen_nacre3_fn
            nacre_name = 'NACRU'
         end select

         beta%ref  = "NETB1"
         oda%ref   = "ODA94"
         lmp%ref   = "LMP00"
         nacre%ref = nacre_name
         illi%ref  = "ILI01"

#ifndef PPN
         if (ipid == 0) then
#endif

            open(unit = netgen_beta_fh  , file = netgen_beta_fn , status = 'old')
            open(unit = netgen_oda_fh   , file = netgen_oda_fn  , status = 'old')
            open(unit = netgen_lmp_fh   , file = netgen_lmp_fn  , status = 'old')
            open(unit = netgen_nacre_fh , file = nacrefn        , status = 'old')
            open(unit = netgen_illi_fh  , file = netgen_illi_fn , status = 'old')

            call read_netgen_file( beta    , netgen_beta_fh  )
            call read_netgen_file( oda     , netgen_oda_fh   )
            call read_netgen_file( lmp     , netgen_lmp_fh   )
            call read_netgen_file( nacre   , netgen_nacre_fh )
            call read_netgen_file( illi    , netgen_illi_fh  )

            close(netgen_beta_fh)
            close(netgen_oda_fh)
            close(netgen_lmp_fh)
            close(netgen_nacre_fh)
            close(netgen_illi_fh)

#ifndef PPN
         end if

         call netgen_broadcasts()
#endif

         allocate ( beta  % rates ( beta  % ireac ) ) 
         allocate ( illi  % rates ( illi  % ireac ) ) 
         allocate ( oda   % rates ( oda   % ireac ) ) 
         allocate ( lmp   % rates ( lmp   % ireac ) ) 
         allocate ( nacre % rates ( nacre % ireac ) ) 

         beta  % rates = ZERO
         illi  % rates = ZERO
         oda   % rates = ZERO
         lmp   % rates = ZERO
         nacre % rates = ZERO
   end subroutine netgen_init


   !> subroutine: netgen_interpolate_rates
   !> 
   !> This routine provides the reaction rates netgenv(i,3) (i=1,nreac)
   !> for temperature and (electron) density
   !> from a linear interpolation of the log of the grid-point rates read
   !> from the Netgen files:
   !>
   !> NACRE
   !> -----
   !>      netgen_nacre_log_100_low.txt,
   !>      netgen_nacre_log_100_rec.txt,
   !>      netgen_nacre_log_100_up.txt.
   !>
   !> are low, rec and upper limit rates of NACRE compilation, in a logaritmic
   !> grid of 100 temperatures between T8=0.01 and T8=100.
   !>
   !> ILIADIS+ 2001
   !> -------------
   !>       netgen_iliadis2001_log_100.txt
   !>
   !> provides the rates of Iliadis et al. 2001 compilation,
   !> in a logaritmic grid of 100 temperatures between T8=0.01 and T8=100.
   !>
   !> BETA DECAYS
   !> -----------
   !>       netgen_betarecom_log_30.txt
   !>
   !> contains rec beta decay rates of NETGEN websource, in a logaritmic grid
   !> of 30 temperatures between T8=0.01 and T8=100. density dependence is
   !> included if required, in a different range according to the reference.
   !>
   !> Note that the factorials accounting for identical particles have
   !> already been included in v(i).
   !> To obtain the evolution dYj/dt of species j, simply multiply v(i) by
   !> the molar fractions of the reacting species j1 + j2: Yj1^n1 Yj2^n2,
   !> where n1, n2 are the stoechiometric factors.
   !>
   !> The number of grid points and reaction rates are determined
   !> automatically, and stored in variables IREAC and KGRID and NDATA
   !>
   !> The rates on grid points are stored in the array
   !>     vgrid(ngrid,nre1,11)
   !> where the last index refers to the density grid for the weak rates
   !> dependent upon electron number density:
   !> vgrid(ngrid,nre,i) with i=1,11
   !> corresponds to the rate at different electron densities N_e as given by
   !> the variable aNegrid(nre,i)
   !> For the other reactions (with flags < 0),
   !> only vgrid(ngrid,nre,1) is relevant
   !>
   !> Coding of reactions:
   !>
   !>   n1 to n4: stoechiometric factors
   !>   z1 to z4: chemical symbol
   !>   a1 to a4: atomic mass
   subroutine netgen_interpolate_rates( t9in, rhoin, yein, ad )
         use utils, only: r8, wallclocktime, fmtprofile
         use constants, only: avogadro
         implicit none

         type(netgen_data) :: ad !< active netgen data

         ! these were set all to 1, I think because they are already included in the netgen rate (?)
         real(r8), parameter :: &
               factorial(0:4) = [ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8 ]

         real(r8), intent(in) :: t9in, rhoin, yein
         ! netgenv is each rate interpolated at the requested T, rho
         real(r8) :: t9, t8, rho, ye, ane_in, ane, l10ratelo, l10ratehi, l10rate, &
               h, b, deltagrid, vbetalow, vbetahigh, tti1, rhop2, rhop3, rhoye, mxr, oldmxr,&
               flg, oldflg, anehi, anelo, dane, lane, lanelo, lanehi, dlane, deltalgrid
         real(r8), allocatable :: netgenv(:)
         integer :: i, j, k, klo, khi, jlo, jhi, ndt, oldndt

         allocate(netgenv(ad%ireac))
         netgenv = ZERO

         ! pass
         t9 = t9in; rho = rhoin; ye = yein

         t8   = t9 * 10._r8 ; rhoye = rho * ye ; ane_in = rhoye * avogadro
         rhop2 = rho * rho ; rhop3 = rho * rho * rho
         mxr = -ONE; ndt = -1; flg = -ONE

         ! ane is the electron density

         ! clip and locate temperature in reaction rate table
         t8 = min(max(t8, ad%tgrid(1)), ad%tgrid(ad%kgrid))
         call bsearch_r8( ad%tgrid(1:ad%kgrid), t8, k)
         klo = k ; khi = klo + 1

         h = ad%l10t8grid(khi) - ad%l10t8grid(klo)
         b = ( log10(t8) - ad%l10t8grid(klo) ) / h

         do i = 1, ad%ireac

            if (ad%vgrid(klo,i,1) <= ZERO .or. ad%vgrid(khi,i,1) <= ZERO) cycle

            if (ad%flag(i) <= ZERO) then
               ! reaction is not electron-density-dependent beta-decay;
               ! linear interpolation of log(rate) in log(t8)
               l10ratelo   = ad%l10vgrid(klo,i,1)
               l10ratehi   = ad%l10vgrid(khi,i,1)
               l10rate = l10ratelo + b * (l10ratehi - l10ratelo)

               select case(nint(ad%flag(i)))
               case(-14)
                  ! electron capture on Be7
                  netgenv(i) = 10._r8 ** l10rate * rhoye
                  if ( t8 < 0.01_r8 ) netgenv(i) = min( netgenv(i), 1.51e-7_r8 / rhoye )
               case(-13)
                  ! electron capture
                  netgenv(i) = 10._r8 ** l10rate * rhoye
               case(-11)
                  ! photodisintegration or beta-decay
                  netgenv(i) = 10._r8 ** l10rate
               case(-10)
                  ! two-particle reaction
                  ! if identical particles: factorials!
                  netgenv(i) = 10._r8 ** l10rate * rho / factorial(ad%n1(i)) / factorial(ad%n2(i))
               case(-100)
                  ! three-particle reactions
                  ! if identical particles: factorials!
                  netgenv(i) = 10._r8 ** l10rate * rhop2 / factorial(ad%n1(i)) / factorial(ad%n2(i))
               case(-200)
                  ! four-particle reactions
                  ! if identical particles: factorials!
                  netgenv(i) = 10._r8 ** l10rate * rhop3 / factorial(ad%n1(i)) / factorial(ad%n2(i))
               case default
                  netgenv(i) = ZERO
               end select

            else
               ! the rate has a tabulated density dependence, and should be interpolated
               ! locate the position of current rho within the density grid aNegrid

               ! TODO: I think this need not be done every single rate, but need to be certain
               ! all rates of a type have the same
               ! aNegrid; probably that is the case and this can be moved outside of the loop
               ! TODO: unify interpolation routines, and replace linear with cubic or quintic
               ! splines

               ! clip and find electron density in rate table if dimensions have changed
               oldflg = flg; flg = ad%flag(i)
               oldndt = ndt; ndt = ad%ndata(i)
               oldmxr = mxr; mxr = ad%aNegrid(i,ndt)
               if (flg /= oldflg .or. ndt /= oldndt .or. oldmxr /= mxr) then
                  !> clip electron density to table
                  ane = min(max(ane_in, flg), mxr)
                  call bsearch_r8( ad%aNegrid(i,1:ndt), ane, j)
                  jlo = j ; jhi = jlo + 1
                  anelo = ad%aNegrid(i,jlo); anehi = ad%aNegrid(i,jhi)
                  dane = ane - anelo
                  deltagrid = anehi - anelo
                  ! log10:
                  lane = log10(ane)
                  lanelo = log10(anelo); lanehi = log10(anehi)
                  dlane = lane - lanelo
                  deltalgrid = lanehi - lanelo
               end if

               ! linear interpolation in log10(temperature) of log10(rate) for two density grid points
               vbetalow  = ad%l10vgrid(klo,i,jlo) + b * (ad%l10vgrid(khi,i,jlo) - ad%l10vgrid(klo,i,jlo))
               vbetahigh = ad%l10vgrid(klo,i,jhi) + b * (ad%l10vgrid(khi,i,jhi) - ad%l10vgrid(klo,i,jhi))

               ! linear interpolation of log10(rate) in electron density
               !netgenv(i) = (vbetahigh - vbetalow) / deltagrid * dane + vbetalow
               !netgenv(i) = 10._r8 ** netgenv(i)

               ! linear interpolation of log10(rate) in log10(electron density)
               !netgenv(i) = (vbetahigh - vbetalow) / deltalgrid * dlane + vbetalow
               !netgenv(i) = 10._r8 ** netgenv(i)

               ! linear interpolation of rate in electron density
               vbetalow  = TEN ** vbetalow
               vbetahigh = TEN ** vbetahigh
               netgenv(i) = (vbetahigh - vbetalow) / deltagrid * dane + vbetalow

            endif
         enddo

         netgenv = max(netgenv,1e-99_r8)

         ! write rates to correct array
         ad%rates = netgenv
         deallocate(netgenv)
   end subroutine netgen_interpolate_rates



   subroutine read_netgen_file( ad, file_handle )
         use reaction_info
         use nuc_data, only: atomic
         implicit none

         type(netgen_data) :: ad

         integer :: file_handle

         real(r8), pointer :: &
               tgrid(:), &                         ! temperature grid
               l10t8grid(:), &                     ! log10 temperature (t8) grid
               aNegrid(:,:), &                     ! electron density grid
               vgrid(:,:,:), &                     ! rates
               l10vgrid(:,:,:), &                  ! log10 rates
               flag(:), &                          ! still do not know
               Qrad(:), &
               Qnu(:)

         ! netgen's own variables used to read the files
         character(len=132) :: longline(7), longlineblank(7)
         integer            :: leng, l, i, j, i1, ifour, ll
         character(len=8)   :: aflag(11)
         integer            :: nn1(1), nn2(1), nn3(1), nn4(1)
         character(len=6)   :: zz1(1), zz2(1), zz3(1), zz4(1)
         real(r8)           :: QQrad(11), QQnu(11)
         character(len=2)   :: netgenel(2)
         integer            :: netgenaa(2)
         character(len=1)   :: dummy

         ! netgen_reac_type: 1=n,g; 2=g,n; 3=n,p; 4=n,a; 6=g,p; 10=g,a 
         integer        :: netgen_reac_type(nre1)

         ! initialise
         ad%z1(:)    = 'XXXXXX' ; ad%z2(:)    = 'XXXXXX' ; ad%z3(:)    = 'XXXXXX' ; ad%z4(:)    = 'XXXXXX'
         ad%n1(:)    = 0        ; ad%n2(:)    = 0        ; ad%n3(:)    = 0        ; ad%n4(:)    = 0
         ad%a1(:)    = 0._r8    ; ad%a2(:)    = 0._r8    ; ad%a3(:)    = 0._r8    ; ad%a4(:)    = 0._r8
         ad%zt(:) = 0           ; ad%zf(:) = 0           ; ad%at(:) = 0           ; ad%af(:) = 0
         ad%elt(:) = 'XX'       ; ad%elf(:) = 'XX'

         ! read the file
         ad%ireac = 0

         do j = 1, 132
            longlineblank(:)(j:j)=' '
         enddo

         do ll = 1, 50000

            ! initialize the string
            longline(:) = longlineblank(:)

            ! read the header
            read(file_handle,1100,end=9999) (longline(j),j=1,7)
            ad%ireac = ad%ireac + 1


            ! l = number of data records on line
            leng = index(longline(5),'         ')
            if(leng == 0) leng = 132
            l = (leng - 11) / 11
            ad%ndata(ad%ireac) = l


            read(longline(7),1200) (aflag(j),j = 1, l)
            ad%aNegrid(ad%ireac,1) = 1._r8


            ! l > 1: flag corresponds to electron densities
            ! read electron densities
            if (l == 1) then
               ad%aNegrid(ad%ireac,1) = 1._r8
            else if (l > 1) then
               read(longline(7),1210) (ad%aNegrid(ad%ireac,j),j=1,l)
               ad%flag(ad%ireac) = ad%aNegrid(ad%ireac,1)
            endif

            read(longline(1),1201,end=1) nn1(1)
            1 read(longline(1),1202,end=2) zz1(1)
            2 read(longline(2),1201,end=3) nn2(1)
            3 read(longline(2),1202,end=4) zz2(1)
            4 read(longline(3),1201,end=5) nn3(1)
            5 read(longline(3),1202,end=6) zz3(1)
            6 read(longline(4),1201,end=7) nn4(1)
            7 read(longline(4),1202,end=8) zz4(1)
            8 read(longline(5),1205,end=11) (QQrad(j),j=1,l)
            11 read(longline(6),1205,err=20,end=12) (QQnu(j), j=1,l)

            12 continue

            ! *** this is done to be able to have neutron and protons as main parent or
            !     dauther specie without having to modify in the file their name.
            if (zz1(1) == 'PROT ') zz1(1) = 'H   1'
            if (zz4(1) == 'PROT ') zz4(1) = 'H   1'
            if (zz1(1) == 'NEUT ') zz1(1) = 'nn  1'
            if (zz4(1) == 'NEUT ') zz4(1) = 'nn  1'

            read(zz1(1),1000,end=14) netgenel(1),netgenaa(1)
            14 read(zz4(1),1000,end=15) netgenel(2),netgenaa(2)

            15 continue

            1000              format(A2,I3,1x)

            ! read the data

            read(file_handle,*)
            read(file_handle,*)

            ! check the grid size
            ! read the first reaction (ad%ireac=1) to find out kgrid
            if (ad%ireac == 1) then

               do i = 1, ngrid
                  read(file_handle,1206,err=21) dummy, ad%tgrid(i), (ad%vgrid(i,ad%ireac,j),j=1,l)

                  if (dummy(1:1)  ==  '#') goto 13
                  ad%kgrid = i

               enddo

               do i = 2,ngrid
                  do j = 1, l
                     if (ad%vgrid(i,ad%ireac,j) == 0) then
                        ad%vgrid(i,ad%ireac,j) = ad%vgrid(i-1,ad%ireac,j)
                     end if
                  end do
               end do

            else

               ! once defined kgrid in the previous step, now read all the
               ! other reactions data.
               ! here we are in the do loop ll and in the do loop nsp (the
               ! ad%ireac one!).

               do i = 1, ad%kgrid
                  read(file_handle,1206,end = 13,err=21) dummy, &
                        ad%tgrid(i),(ad%vgrid(i,ad%ireac,j),j = 1,l)
               enddo

               do i = 2, ad%kgrid
                  do j = 1, l
                     if (ad%vgrid(i,ad%ireac,j) == 0) then
                        ad%vgrid(i,ad%ireac,j) = ad%vgrid(i-1,ad%ireac,j)
                     end if
                  end do
               end do

            endif

            read(file_handle,*)

            13 continue

            ad%n1(ad%ireac)      = nn1(1)
            ad%n2(ad%ireac)      = nn2(1)
            ad%n3(ad%ireac)      = nn3(1)
            ad%n4(ad%ireac)      = nn4(1)
            ad%z1(ad%ireac)      = zz1(1)
            ad%z2(ad%ireac)      = zz2(1)
            ad%z3(ad%ireac)      = zz3(1)
            ad%z4(ad%ireac)      = zz4(1)
            ad%Qrad(ad%ireac)    = int(QQrad(1))
            ad%Qnu(ad%ireac)     = int(QQnu(1))

            ad%elt(ad%ireac)  = netgenel(1)
            ad%elf(ad%ireac)  = netgenel(2)
            ad%at(ad%ireac)   = netgenaa(1)
            ad%af(ad%ireac)   = netgenaa(2)

            ! z are calcultated istead of element label: useful for merging in general network
            call atomic(nre1,ad%elt,ad%elf,ad%zt,ad%zf,ad%ireac)

            ! I try to define here ntransv, once I have also atomic number.
            ! not included: 1=n,g; 2=g,n; 3=n,p; 4=n,a; 6=g,p; 10=g,a
            ! maybe I do not need ntransv!! let's see....
            if(zz2(1)(1:4).eq.'PROT'.and.zz3(1)(1:5).eq.'OOOOO')then
               netgen_reac_type(ad%ireac) = i_pg
            else if(zz2(1)(1:4).eq.'PROT'.and.zz3(1)(1:4).eq.'NEUT')then
               netgen_reac_type(ad%ireac) = i_pn
            else if(zz2(1)(1:4).eq.'PROT'.and.zz3(1)(1:5).eq.'HE  4')then
               netgen_reac_type(ad%ireac) = i_pa
            else if(zz2(1)(1:5).eq.'HE  4'.and.zz3(1)(1:5).eq.'OOOOO')then
               netgen_reac_type(ad%ireac) = i_ag
            else if(zz2(1)(1:5).eq.'HE  4'.and.zz3(1)(1:4).eq.'NEUT')then
               netgen_reac_type(ad%ireac) = i_an
            else if(zz2(1)(1:5).eq.'HE  4'.and.zz3(1)(1:4).eq.'PROT')then
               netgen_reac_type(ad%ireac) = i_ap
            else if(zz2(1)(1:5).eq.'OOOOO'.and.zz3(1)(1:5).eq.'OOOOO')then
               if(ad%zt(ad%ireac).lt.ad%zf(ad%ireac))then
                  netgen_reac_type(ad%ireac) = i_bm
               else
                  netgen_reac_type(ad%ireac) = i_ec
               end if
            else
               ! case for nacre, for initial reactions different from the ones
               ! listed above...(case 4 does not correpond to (n,a) here!!!!!)
               ! THEN WHY IS IT SET TO 4?!?!?!?!?!?!?!?!!>!>!>!>!?>! ARGHHHHH
               netgen_reac_type(ad%ireac) = 4
            end if
            ! ***
            if(aflag(1)(5:8).eq.'----') then
               ad%flag(ad%ireac) = -200.
            else if(aflag(1)(6:8).eq.'---') then
               ad%flag(ad%ireac) = -100.
            else if(aflag(1)(7:8).eq.'--') then
               ad%flag(ad%ireac) = -10.
            else if(aflag(1)(6:8).eq.'+++') then
               ad%flag(ad%ireac) = -14.
            else if(aflag(1)(7:8).eq.'++') then
               ad%flag(ad%ireac) = -13.
            else if(aflag(1)(8:8).eq.'+') then
               ad%flag(ad%ireac) = -11.
            endif


            if(zz1(1)(1:4).eq.'NEUT'.or.zz1(1)(1:4).eq.'PROT') then
               ad%a1(ad%ireac) = ONE
            else if(zz1(1)(1:5).eq.'OOOOO') then
               ad%a1(ad%ireac) = ZERO
            else if(zz1(1)(1:4).eq.'DEUT') then
               ad%a1(ad%ireac) = TWO
            else if(zz1(1)(1:4).eq.'TRIT') then
               ad%a1(ad%ireac) = THREE
            else
               read(zz1(1)(3:5),'(i3)') i1
               ad%a1(ad%ireac) = float(i1)
            endif

            if(zz2(1)(1:4).eq.'NEUT'.or.zz2(1)(1:4).eq.'PROT') then
               ad%a2(ad%ireac) = ONE
            else if(zz2(1)(1:5).eq.'OOOOO') then
               ad%a2(ad%ireac) = ZERO
            else if(zz2(1)(1:4).eq.'DEUT') then
               ad%a2(ad%ireac) = TWO
            else if(zz2(1)(1:4).eq.'TRIT') then
               ad%a2(ad%ireac) = THREE
            else if(zz2(1)(1:5).eq.'HE  4') then
               ad%a2(ad%ireac) = FOUR
            else if(zz2(1)(1:5).eq.'HE  3') then
               ad%a2(ad%ireac) = THREE
            endif
            if(zz3(1)(1:4).eq.'NEUT'.or.zz3(1)(1:4).eq.'PROT') then
               ad%a3(ad%ireac) = ONE
            else if(zz3(1)(1:5).eq.'OOOOO') then
               ad%a3(ad%ireac) = ZERO
            else if(zz3(1)(1:4).eq.'DEUT') then
               ad%a3(ad%ireac) = TWO
            else if(zz3(1)(1:4).eq.'TRIT') then
               ad%a3(ad%ireac) = THREE
            else if(zz3(1)(1:5).eq.'HE  4') then
               ad%a3(ad%ireac) = FOUR
            else if(zz3(1)(1:5).eq.'HE  3') then
               ad%a3(ad%ireac) = THREE
            endif
            if(zz4(1)(1:4).eq.'NEUT'.or.zz4(1)(1:4).eq.'PROT') then
               ad%a4(ad%ireac) = ONE
            else if(zz4(1)(1:5).eq.'OOOOO') then
               ad%a4(ad%ireac) = ZERO
            else if(zz4(1)(1:4).eq.'DEUT') then
               ad%a4(ad%ireac) = TWO
            else if(zz4(1)(1:4).eq.'TRIT') then
               ad%a4(ad%ireac) = THREE
            else
               read(zz4(1)(3:5),'(i3)') ifour
               ad%a4(ad%ireac) = real(ifour,r8)
            endif

            ! ^_^ The writing of ntrans at the moment, unforunately, has to be done so, because of the
            ! complicated indexing

            select case(file_handle)
            case(netgen_beta_fh)
               netgen_beta_ntrans(  ad%at(ad%ireac), ad%zt(ad%ireac),1, netgen_reac_type(ad%ireac)) = ad%ireac
            case(netgen_oda_fh)
               netgen_oda_ntrans(   ad%at(ad%ireac), ad%zt(ad%ireac),1, netgen_reac_type(ad%ireac)) = ad%ireac
            case(netgen_lmp_fh)
               netgen_lmp_ntrans(   ad%at(ad%ireac), ad%zt(ad%ireac),1, netgen_reac_type(ad%ireac)) = ad%ireac
            case(netgen_illi_fh)
               netgen_illi_ntrans(  ad%at(ad%ireac), ad%zt(ad%ireac),1, netgen_reac_type(ad%ireac)) = ad%ireac
            case(netgen_nacre_fh)
               netgen_nacre_ntrans( ad%at(ad%ireac), ad%zt(ad%ireac),1, netgen_reac_type(ad%ireac)) = ad%ireac
            end select

         enddo

         1100    format(4(a132,//),3(a132,/))
         1200    format(14x,11(a8,3x))
         1201    format(14x,11(i1,10x))
         1202    format(14x,11(2x,a6,3x))
         ! 1204    format(i1,1x,a6,'( ',i1,1x,a5,',',1x,i1,1x,a5,')',2x,i1,1x,a6)
         1205    format(14x,11(f8.3,3x))
         1206    format(a1,2x,f8.4,11(1x,e10.4))
         1210    format(14x,11(e8.1,3x))

         9999        continue

         if (ad%tgrid(1) > 0.01d0) print*,'warning, enlarge your netgen table' // &
              'with file handle =', file_handle, ' to lower temperatures: tgrid(1) =', ad%tgrid(1)
         if (ad%tgrid(ad%kgrid) < 100.d0) print*,'warning, enlarge your netgen table' // &
               'with file handle =', file_handle, ' to higher temperatures'

         ! store log10 of all rates, because the linear interpolation is done in log, and we can
         ! avoid taking expensive logarithms multiple times this way

         ad%l10vgrid(:,:,:) = -1.e99_r8
         ad%l10t8grid(:) = -1.e99_r8
         where ( ad%vgrid > 0._r8 )
            ad%l10vgrid(:,:,:) = log10(ad%vgrid(:,:,:))
         end where
         where ( ad%tgrid > 0._r8 )
            ad%l10t8grid(:)     = log10(ad%tgrid(:))
         endwhere

         return

         ! ***  branch here if NaN was detected in Qnu line
         20 continue

         write(*,7000) longline(1)(15:21),longline(2)(15:21), &
               longline(3)(15:21),longline(4)(15:21)

         7000    format('** WARNING: NaN have been detected among Qnu values', &
               ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/, &
               ' This may be because Qnu is density-dependent.', &
               ' The Netgen log file provides a complete table', &
               ' of Qnu(T,rho).',/,' Edit the vit.dat file to remove', &
               ' NaNs and try again')

         stop

         !        branch here if NaN was detected among V
         21 continue

         write(*,7001) longline(1)(15:21),longline(2)(15:21), &
               longline(3)(15:21),longline(4)(15:21)

         7001    format('** WARNING: NaN has been detected in reaction rate', &
               ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/, &
               'This may be because your grid extends beyond the', &
               ' valid data range for this reaction', &
               ' (see Netgen log file)',/, &
               ' Either change your temperature or density grid,', &
               ' or remove NaN in vit.dat by extrapolating the', &
               ' existing data (risky!)')

         stop

   end subroutine read_netgen_file


   integer(i4) function netgen_rtype(ad, i)
         use reaction_info
         implicit none
         type(netgen_data) :: ad
         character(6) :: z2, z3
         integer(i4) :: i, zt, zf

         z2 = ad % z2(i); z3 = ad % z3(i); zt = ad % zt(i); zf = ad % zf(i)
         netgen_rtype = -1

         if (z2 == 'OOOOO ' .and. z3 == 'OOOOO ') then
               if (zt < zf) then
                  ! beta decay
                  netgen_rtype = i_bm
               else if (zt > zf) then
                  ! electron capture
                  netgen_rtype = i_ec
               end if
         else if (z2(1:4) == 'PROT' .and. z3(1:5) == 'OOOOO') then
            if (zt < zf) then
               netgen_rtype = i_pg
            end if
         else if (z2(1:4) == 'PROT' .and. z3(1:4) == 'NEUT') then
            if (zt < zf) then
               netgen_rtype = i_pn
            end if
         else if (z2(1:4) == 'PROT' .and. z3(1:5) == 'HE  4') then
            if (zt > zf) then
               netgen_rtype = i_pa
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:5) == 'OOOOO') then
            if (zt < zf) then
               netgen_rtype = i_ag
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:4) == 'NEUT') then
            if (zt < zf) then
               netgen_rtype = i_an
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:4) == 'PROT') then
            if (zt < zf) then
               netgen_rtype = i_ap
            end if
         else if (z2(1:4) == 'PROT' .and. z3(1:5) == 'OOOOO') then
            if (zt < zf) then
               netgen_rtype = i_pg
            end if
         else if (z2(1:4) == 'PROT' .and. z3(1:4) == 'NEUT') then
            if (zt < zf) then
               netgen_rtype = i_pn
            end if
         else if (z2(1:4) == 'PROT' .and. z3(1:5) == 'HE  4') then
            if (zt > zf) then
               netgen_rtype = i_pa
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:5) == 'OOOOO') then
            if (zt < zf) then
               netgen_rtype = i_ag
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:4) == 'NEUT') then
            if (zt < zf) then
               netgen_rtype = i_an
            end if
         else if (z2(1:5) == 'HE  4' .and. z3(1:4) == 'PROT') then
            if (zt < zf) then
               netgen_rtype = i_ap
            end if
         end if
   end function netgen_rtype


#ifndef PPN
   subroutine netgen_broadcasts()

         call broadcast_netgen_data(beta)
         call broadcast_netgen_data(oda)
         call broadcast_netgen_data(lmp)
         call broadcast_netgen_data(illi)
         call broadcast_netgen_data(nacre)
         call broadcast(netgen_beta_ntrans)
         call broadcast(netgen_nacre_ntrans)
         call broadcast(netgen_illi_ntrans)
         call broadcast(netgen_oda_ntrans)
         call broadcast(netgen_lmp_ntrans)

   end subroutine netgen_broadcasts
#endif



end module netgen
