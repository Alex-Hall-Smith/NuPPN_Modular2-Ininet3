module vital

   use utils, only: r8, table_input
   use constants
   use array_sizes
   use rates, only : v, dvdt, dvdd, k1, k2, k3, k4, k5, k6, k7, k8
   use frame_knobs
   use nuc_data, only: ispe, considerisotope, considerreaction, t_half, t_half_unit, zis
   use physics_knobs, only: tbetamin, getderivs

   implicit none

   ! ^_^ see end of module for list of VITAL reactions

   ! ^_^ BLOCK DATA for o18(a,g) reaction
   integer, parameter :: NG = 53
   real(r8) :: TG(NG),O18AG(NG)
   
   ! ^_^ rate tables moved to init for reading
   !     new Ne22+a by Michael Wiescher, Joachim Gorres, G. Imbriani, J.
   !     deBoer and Mary Beard (Spring 2014)
   !     First implemented by Nobuya (Nishimura N. et al. 2014 AIPC)

   integer, parameter:: ndat_ne22 = 51

   integer:: ne22an, ne22ag, idt, i_ne22,i_c12
   integer, parameter :: ndat_c12ag = 39
   real(r8) :: &
         t9_in(1:ndat_ne22), &
         rt_ne22an(1:ndat_ne22), &
         rt_ne22ag(1:ndat_ne22), &
         rt_c12ag(1:ndat_c12ag), &
         rt_ne22an_in(1:ndat_ne22,1:3), &
         rt_ne22ag_in(1:ndat_ne22,1:3), &
         rt_c12ag_in(1:ndat_c12ag,1:3)
   real(r8) :: h_intpl
   real(r8) :: rate_ne22an_mw,rate_ne22ag_mw,rate_ne22an_lo, &
         rate_ne22ag_lo,rate_ne22an_nd,rate_ne22ag_nd, &
         rate_c12ag_nd
   logical:: ne22_michael = .false.
   logical:: ne22_longland = .false.
   logical:: ne22_nd2015 = .false.
   logical:: c12ag_jdb2016 = .false.

   integer, parameter :: nc12t = 50, nc12 = 3
   real(r8) :: t9c12(nc12t), c12tab(nc12,nc12t)

   ! ^_^ to be burned violently:
   common/mwc12c12/t9c12,c12tab

   ! formats
   character(len=256), parameter :: &
         fmtspecies  = "(I4,1X,A5,1X,2(F4.0),1X,L1,1x,1ES10.3,1x,A1)", &
         fmtreaction = "(I4,1X,L1,1X,I2,1X,A5,2X,I2,1X,A5,1X,I2,1X,A5,2X,I2,1X,A5,"//&
                       "3X,F7.3,3x,A5,2x,A5,2x,I2,2x,1ES10.3)"

   private

   public vital_init, vital_rates_derivs, vital_index, read_physics_input_data


contains



   subroutine read_physics_input_data(an, zn, qi, nrcp, nvcp)
         use communication
         use reaction_info, only: lab, labb, ilabb, rfac

         integer :: i, l, nnnz, nnnr, ispec(4), ndum, nrcp, nvcp
         character(len=80) :: cdummy
         character(len=5) :: spe1, spe2, spe3, spe4
         real(r8) :: an(nsp), zn(nsp), qi(nre)

         ! read VITAL species

         if (master) then
            read (2, '(a80)') cdummy
            do i = 1, nvcp
               read(2, fmtspecies) nnnz, zis(i), an(i), zn(i), considerisotope(i), t_half(i), t_half_unit(i)
            end do
         end if

#ifndef PPN
         call broadcast(an); call broadcast(zn); call broadcast(t_half)
         call broadcast(considerisotope)
         call broadcast_ch_arr(zis); call broadcast_ch_arr(t_half_unit)
#endif

         ! compare half lives for VITAL species to tbetamin parameter

         do i = 1, nvcp
            if (.not. considerisotope(i)) cycle

            ! convert half life to seconds
            select case(t_half_unit(i))
            case('m')
               t_half(i) = t_half(i) * min2sec
            case('h')
               t_half(i) = t_half(i) * hrs2sec
            case('d')
               t_half(i) = t_half(i) * day2sec
            case('y')
               t_half(i) = t_half(i) * yrs2sec
            end select

            if (master) then
               if (t_half(i) < tbetamin) then
                  write(*,*) 'rnetw2007: in VITAL network 1 half life < tbetamin'
                  write(*,*) 'isotope=',zis(i),'half life=',t_half(i)
                  write(*,*) 'terrestrial decay rate activated in vital'
               end if
            end if
         end do

         ! read VITAL reactions

         if (master) then
            ndum = 0
            do l = 1, nrcp
               read(2, fmtreaction) nnnr, considerreaction(l), &
                     k2(l), spe1, k4(l), spe2, k6(l), spe3, k8(l), spe4, &
                     qi(l), lab(l), labb(l), ilabb(l), rfac(l)

               k1(l) = ispe(spe1)
               k3(l) = ispe(spe2)
               k5(l) = ispe(spe3)
               k7(l) = ispe(spe4)

               ! make sure all reactants of reactions have been selected as isotopes
               if (considerreaction(l)) then
                  ispec(:) = [k1(l), k3(l), k5(l), k7(l)]
                  do i = 1, 4
                     if (.not. considerisotope(ispec(i)) .and. ispec(i) /= ispe("OOOOO")) then
                        ndum = ndum + 1
                        write(*,*)'ERROR in rnetw2007: '
                        write(*,*) 'isotope ', zis(ispec(i)), ' occurs in reaction ', l, &
                              'but is not selected.'
                     endif
                  end do
               else
                  qi(l) = ZERO
               end if

               if (ndum > 0) then
                  write(*,*) 'BAD WARNING in rnetw2007:'
                  write(*,*) 'reactions are selected (T) for which not all species are selected (T)'
               end if
            end do
            close(2)
         end if

#ifndef PPN
         call broadcast(k1)         ; call broadcast(k2)          ; call broadcast(k3)
         call broadcast(k4)         ; call broadcast(k5)          ; call broadcast(k6)
         call broadcast(k7)         ; call broadcast(k8)          ; call broadcast(considerreaction)
         call broadcast(qi)         ; call broadcast(ilabb)       ; call broadcast(rfac)
         call broadcast_ch_arr(lab) ; call broadcast_ch_arr(labb)
#endif
   end subroutine read_physics_input_data



   subroutine vital_init()

         integer  :: i
         real(r8) :: temp

         ! ^_^ read in deBoer (2016) rate, which was previously being done EVERY
         !     TIMESTEP FROM A FILE!!!!!???

         open(table_input, file = '../NPDATA/c12ag_jdb16.dat', action = 'read')

         read(table_input,*)
         read(table_input,*)

         do i = 1, ndat_c12ag-1
            read(table_input,*) t9_in(i), rt_c12ag_in(i,1:3)
         end do
         close(table_input)

         ! ^_^ read in c12c12 (Wiescher??) rate, which was *also* being read in
         !     EVERY TIMESTEP FROM A FILE!!!?????

         open(table_input, file = "../NPDATA/12C+12Crate_new.tex", status = 'old')

         do i = 1, 2
            read(table_input, *)
         end do

         do i = 1, nc12t
            read(table_input, '(10(ES9.2, 3x))') t9c12(i), temp, temp, &
                  c12tab(1,i), temp, temp, c12tab(2,i), temp, temp, c12tab(3,i)
         end do
         close(table_input)

         ! some legacy block data
         TG(:) = [ 0.01D0,0.02D0,0.03D0,0.04D0,0.05D0,0.06D0,0.07D0, &
               0.08D0,0.09D0,0.10D0,0.11D0,0.12D0,0.13D0,0.14D0, &
               0.15D0,0.16D0,0.17D0,0.18D0,0.19D0,0.20D0,0.21D0, &
               0.22D0,0.23D0,0.24D0,0.25D0,0.26D0,0.27D0,0.28D0, &
               0.29D0,0.30D0,0.31D0,0.32D0,0.33D0,0.34D0,0.35D0, &
               0.36D0,0.37D0,0.38D0,0.39D0,0.40D0,0.41D0,0.42D0, &
               0.43D0,0.44D0,0.45D0,0.46D0,0.47D0,0.48D0,0.49D0, &
               0.50D0,0.60D0,0.80D0,1.00D0 ]

         O18AG(:) = [ 1.22D-62,4.04D-51,4.43D-41,8.62D-34,1.89D-29, &
               1.40D-26,1.52D-24,5.02D-23,8.86D-22,2.46D-20, &
               9.58D-19,2.41D-17,3.72D-16,3.88D-15,2.94D-14, &
               1.72D-13,8.16D-13,3.25D-12,1.12D-11,3.42D-11, &
               9.45D-11,2.40D-10,5.68D-10,1.27D-09,2.70D-09, &
               5.49D-09,1.08D-08,2.05D-08,3.78D-08,6.78D-08, &
               1.19D-07,2.02D-07,3.37D-07,5.50D-07,8.77D-07, &
               1.37D-06,2.10D-06,3.15D-06,4.64D-06,6.73D-06, &
               9.60D-06,1.35D-05,1.87D-05,2.55D-05,3.43D-05, &
               4.57D-05,6.01D-05,7.82D-05,1.01D-04,1.28D-04, &
               9.31D-04,1.07D-02,4.50D-02 ]

   end subroutine vital_init


   subroutine vital_rates_derivs( rho, t9_in, xii )
         real(r8), intent(in) :: rho, t9_in, xii(nsp)
         real(r8) :: t9
         integer :: nrpcold
         ! *** nrpc of vital
         integer nrcpold
         common/old/nrcpold

         ! clip temperature
         t9 = min(t9_in, 10._r8)

         ! ^_^ rates
         call vital_calculate_rates( rho, t9, xii )

   end subroutine vital_rates_derivs


   subroutine vital_calculate_rates( rho, t9, xii )

         real(r8), intent(in) :: xii(nsp)
         real(r8) :: x(nsp)
         real(r8) :: BRCCN,BRCCP,BRCCA,BRCON,BRCOP,BRCOA,BROON,BROOP,BROOA,T9CN,T9PPMAX,T9,RHO
         real(r8) :: TP13, TP14, TP12, TP15, TP23, TP32, TP38, TP27, TP43, lnt9, t9inv, l10t9
         real(r8) :: TP53,TP35,TP54,TP65,TP72,TM13,TM23,TM45,ALPHA
         real(r8) :: T912,T913,T914,T915,T923,T927,T932,T935,T943,T953,T954,T972,TM12
         real(r8) :: TM32,T9B,T9B13,T9B23,T9B56,T9H,T9H13,T9H23,T9H56,T9X,T9X13,T9X23
         real(r8) :: T9X56,GT92,GT94,GT96,GT98,T9A,TT13,TT56,GT9,TO, &
               GO,TO13,TO56,VALAMG,VMGASI,VNA23PG,VALPSI,VNA23PA,C12C12, &
               FPT9A,FT9A,VO18,V1,V2,F14,F17,F27A,F27B,ARG1,ARG2,G4,FNA23P,FNA23, &
               FNE22A,FNE22B,FNE21,FC,FA,G2,G1,F1,F2,REV,T9A13,T13,T56, &
               TT,VMAX,FABE7,PROTOFAK,T9A56,O16O16,C12O16,T9CO, &
               T9OO,HE4HE4,HE4BE8
         real(r8) :: ROOTT8, a0r,a1r,a2r,a3r,a4r,a5r,a6r
         integer INETW,ICYCL,IPPIV,I,IPROT
         real(r8) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, &
               a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11, &
               i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12, &
               Vn77,Vr177,Vr277,V_n111,V_r111,V_n115,V_r115,&
               V_n116,V_r116
         real(r8) :: r_o16ag
         integer Buch,kunz,imb,deboer
         ! *** for new c12c12
         real(r8) :: xfind,mask,xfac,c12c12w(3),temp
         integer ifindpos(1)

         ! *** nrpc of vital
         integer nrcpold
         common/old/nrcpold
         ! *** NN ------------------------------------------------------------------

         !
         !.. NUCLEAR REACTION RATES
         !.. ( RHO*NAVOGADRO* (SIGMA*VIT)MAXWELL)  ( SEC-1)
         !
         DATA BRCCP,BRCCA/0.35D0,0.65D0/
         DATA BROON,BROOP,BROOA/0.075d0,0.585D0,0.34D0/
         DATA BRCON,BRCOP,BRCOA/0.1D0,0.5D0,0.4D0/
         DATA T9CN/0.30D0/
         DATA T9CO/0.5d0/
         DATA T9OO/0.35d0/
         DATA INETW/1/           ! ALLE RATEN WERDEN AUSGERECHNET,
         ! INETW=0 -> NUR DIE FUER ENERG
         DATA ICYCL/1/           ! ZUSCHALTEN DER MgAl und NeMa Raten
         DATA IPPIV/0/           ! Raten fuer hot H-deficient He3-burning
         DATA T9PPMAX/10._r8/      ! maximale Temperatur fuer PP-Raten
         c12o16   = ZERO
         o16o16   = ZERO
         r_o16ag  = ZERO


         ROOTT8   = (T9 * 10._r8) ** HALF
         TP13     = T9 ** THIRD
         TP14     = T9 ** FOURTH
         TP12     = TP14 * TP14
         TP15     = T9 ** FIFTH
         TP23     = TP13 * TP13
         TP32     = T9 * TP12
         TP38     = T9 ** THREE_EIGTHS
         TP27     = T9 ** TWO_SEVENTHS
         TP43     = T9 * TP13
         TP53     = T9 * TP23
         TP35     = T9 ** THREE_FIFTHS
         TP54     = T9 * TP14
         TP65     = T9 ** SIX_FIFTHS
         TP72     = (T9 * T9 * T9) * TP12
         TM13     = ONE / TP13
         TM12     = ONE / TP12
         TM23     = ONE / TP23
         TM32     = ONE / TP32
         TM45     = ONE / (T9 ** FOUR_FIFTHS)
         T913     = TP13
         T923     = TP23
         T932     = TP32
         lnt9     = log(t9)
         l10t9    = log10(t9)
         t9inv    = ONE / t9

         x(:)     = xii(:)

         ! set all vital rates to zero
         v   (1:nrcpold) = ZERO

         ! -------------------------------------------------
         !               HYDROGEN BURNING
         ! -------------------------------------------------

         IPROT    = ISPE('PROT ')
         PROTOFAK = ONE
         if (.true.) then
            !---------------------------------------
            !      PP-CHAIN

            !---------------------------------------
            !     TEMPERATURE DEPENDENT DECAYS
            !---------------------------------------

            !- BE7 (E-,NU+G) LI7................................. CF88

            FABE7 = HALF * (ONE + XII(IPROT))
            v(5)  = 1.34e-10_r8 / TP12 * (ONE - 0.537_r8 * TP13 + 3.86_r8 * TP23 + &
                  0.0027_r8 * t9inv * EXV(2.515e-03_r8 *t9inv)) * RHO * FABE7
            if (T9 < 1.e-03_r8 .or. T9 > T9PPMAX) v(5) = min(v(5), 1.51e-07_r8)

            !---------------------------------------
            !     TEMPERATURE INDEPENDENT DECAYS
            !---------------------------------------
            !- N13 (E+ NU)  C13
            !            V(91)=1.16D-3
            !C            V(91)=1.9D-11 !1000yt
            !C104   -   1  C  14  +  1  E-    YPS  1  N  14  +  0  NU          ! Karlsruher NC
            !            V(104)=3.84D-12 !ABV,p26
            !C105   -   1  NI 59  +  1  E+      1  CO 59  +  0  NU          ! Karlsruher NC
            !            V(105)=2.93E-13
            !C 115 T  1 NI 63 ( 1 OOOOO, 1 CU 63)  0 OOOOO +   0.0
            !            V(115)=2.200D-10

            !---------------------------------------
            ! *** CF88 and nacre say that below t9 = 0.001 (T = 1.d+6 K)
            ! *** formula are not correct.
            ! *** So, except for Be7 e+ (see v(5) above), for t9 < 0.001
            ! *** I set all the rates of vital = 0.

            if (t9 >= 0.001_r8) then

               IF (T9 > 1.2e-3_r8 .AND. T9 < T9PPMAX) THEN
                  FABE7 = HALF * (ONE + XII(IPROT))
                  !---------------------------------------------------
                  !                 PP I
                  !---------------------------------------------------

                  !-  H(P,E+NU) D                         FILIPPONE 1985
                  !      V(1)=(3.82d-15/TP23*EXV(-3.38/TP13)*(ONE+0.123*TP13+1.09*TP23+
                  !     &     0.938*T9))*RHO*PROTOFAK


                  !-  H(P,E+NU) D                         CF 88

                  V(1) = (4.01e-15_r8 * TM23 * EXV(-3.380_r8 *TM13)* &
                        (ONE + 0.123_r8 * TP13 + 1.09_r8 * TP23 + 0.938_r8 * T9)) &
                        * RHO * PROTOFAK

                  !-  H(P,E+NU) D        (Q=1.442 MeV)         NACRE
                  ! *** better to not use nacre fit, but tables
                  ! *** this seems fine.
                  !      v(1)= (4.08d-15/TP23 * EXV(-3.381d0/TP13)*
                  !     1 (1.d0 + 3.82d0*T9 + 1.51d0*T9**2. + 0.144d0*T9**3.
                  !     2 - 1.14d-2*T9**4.))
                  !     3 * RHO * PROTOFAK

                  !      print *,'H(P,E+NU)D NACRE',v(1)/rho/PROTOFAK

                  !****************************************************************
                  !-  D(P,G) He3...................CF 88
                  V(2) = (2.24D+03 * TM23 * EXV(-3.720_r8 * TM13) * &
                        (ONE + 0.112_r8 * TP13 + 3.38_r8 * TP23 + 2.65_r8 * T9)) * RHO * PROTOFAK

                  !-  D(P,G) He3......(Q=5.493 MeV).......NACRE
                  ! *** better to not use nacre fit, but tables
                  ! *** this seems fine.
                  !      if(t9.le.0.11)then
                  !       V(2)= (1.81D+03/TP23 * EXV(-3.721/TP13) *
                  !     1  (ONE+ 14.3d0*T9 - 90.5d0*T9**2. + 395.d0*T9**3.))
                  !     2  * RHO * PROTOFAK
                  !      else if(t9.gt.0.11)then
                  !       V(2)= (2.58D+03/TP23 * EXV(-3.721/TP13) *
                  !     1  (1. + 3.96d0*T9 + 0.116d0*T9**2.))
                  !     2  * RHO * PROTOFAK
                  !      end if

                  ! *** reverse rate: He3(g,p)D     to be tested!
                  !     V = V(2)/(RHO*PROTOFAK)
                  !     * 1.630d+10*T9**(3./2.)exv(-63.749d0/T9)


                  !****************************************************************

                  !   -HE3 (HE3,2P) HE4 ...........CF 88
                  V(3) = (6.04D+10 * TM23 * EXV(-12.276_r8 * TM13) * (1._r8 + &
                        0.034_r8 * TP13 - 0.522_r8 * TP23 - 0.124_r8 * T9 + 0.353_r8 * &
                        TP43 + 0.213_r8 * TP53)) * RHO

                  !   -HE3 (HE3,2P) HE4 ..(Q=12.859 MeV)......NACRE
                  ! *** better to not use nacre fit, but tables
                  ! *** this seems fine.
                  !      V(3) = (5.59D+10/TP23 * EXV(-12.277/TP13) *
                  !     1 (1. - 0.135d0*T9 + 2.54d-2*T9**2. - 1.29d-3*T9**3.))
                  !     2 * RHO

                  ! *** reverse rate: He4(2p,He3)He3     to be tested!
                  !     V = V(3)/(RHO)
                  !     * 3.392d-10*T9**(-3./2.)exv(-149.23d0/T9)



                  !****************************************************************

                  !---------------------------------------------
                  !                PP II
                  !---------------------------------------------
                  !- HE4 (HE3,G) BE7...............CF 88
                  TT = T9 / (ONE + 4.95e-02_r8*T9)
                  T56 = TT ** FIVE_SIXTHS
                  T13 = TT ** THIRD
                  V(4) = (5.61e+06_r8*T56/TP32 * EXV(-12.826_r8/T13)) * RHO

                  ! *** NACRE formula has a problem at high temperatures.
                  ! *** T9=9 factor 2.36 lower than tables.
                  !- HE4 (HE3,G) BE7.......(Q=1.587 MeV).....NACRE
                  !      v(4) = 5.46d+6/T9**(2./3.) * exv(-12.827d0/T13) *
                  !     1 (1.d0 - 0.307d0*T9 + 8.81d-2*T9**2. -
                  !     2 1.06d-2*T9**3. + 4.46d-4*T9**4.)
                  !     3 * RHO

                  ! *** reverse rate: BE7(g,He3)He4     to be tested!
                  !     V = V(4)/(RHO)
                  !     * 1.113d+10*T9**(3./2.)exv(-18.412d0/T9)

                  !****************************************************************


                  !- LI7 (P,A) HE4                                          FCZ1975
                  !**    V(6)=RHO*(8.04d+08/TP23*EXV(-8.471/TP13-(T9/30.068)**2)*
                  !**   1 (ONE+0.049*TP13+0.230*TP23+0.079*T9-0.027*TP43-0.023*TP53)
                  !**   2 + 1.54d+06/TP32*EXV(-4.479_r8/T9) +1.07d+10/TP32*EXV(-30.443_r8/T9))

                  !- LI7 (P,A) HE4..................... CF 88
                  T9A = T9 / (ONE + 0.759_r8*T9)
                  T9A13 = T9A**THIRD
                  T9A56 = T9A**FIVE_SIXTHS
                  V(6)=(1.096D+09/TP23*EXV(-8.472_r8/TP13) - &
                        4.83D+08*T9A56/TP32*EXV(-8.472_r8/T9A13) + &
                        1.06D+10/TP32*EXV(-30.442_r8/T9))* RHO * PROTOFAK
                  !-----------------------------------------------------------
                  !               PP III
                  !-----------------------------------------------------------
                  !-BE7 (P,G) B8               CF 88
                  V(7) = (3.11e5_r8 / TP23 * EXV(-10.262_r8/TP13) + 2.53e3_r8 / TP32 * EXV(-7.306_r8/T9))* RHO * PROTOFAK

                  !-B8 (G,P) BE7
                  REV = 1.30d+10 * TP32 * EXV(-1.595_r8 / T9) / (RHO)
                  V(9) = REV * V(7)

                  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                  !-B8 (D+NU) 2 ALPHA         ?? is it ok still?                       W1969
                  V(8) = 8.96e-1_r8

                  !-----------------------------------------------------------
                  !               PP IV (HOT H-DEFICIENT HE3-BURNING)
                  !-----------------------------------------------------------
                  IF (IPPIV .EQ. 1) THEN
                     !-BE7 (ALPHA,G) C11                                 NACRE (Angulo etal 1999)
                     ! *** better to not use nacre fit, but tables
                     ! *** included by Falk. To switch to table.
                     IF (T9.LT.2.D0) THEN
                        V(10)=1.29D+10*TM23*EXV(-23.214_r8*TM13-(T9/0.8_r8)**2) &
                              *(ONE -6.47_r8*T9 +19.5_r8*T9*T9 -19.3_r8*T9*T9*T9) &
                              +1.25D+04*TM32*EXV(-6.498_r8/T9) &
                              +1.44D+05*TM32*EXV(-10.177_r8/T9) &
                              +1.63D+04*(T9**0.178_r8)*EXV(-15.281_r8/T9)*RHO
                     ELSE
                        STOP 'BE7 (ALPHA,G) C11: T9 > 2'
                     ENDIF
                     !-LI7 (ALPHA,G) B11                                 NACRE (Angulo etal 1999)
                     ! *** better to not use nacre fit, but tables
                     ! *** included by Falk. To switch to table.
                     V(11)=9.72D+07*TM23 *EXV(-19.163_r8*TM13-(T9/0.4_r8)**2) &
                           *(1.D0+2.84_r8*T9 -7.89_r8*T9*T9) +3.355D+02*TM32 *EXV(-2.959_r8/T9) &
                           +1.04D4*T9**(-0.023_r8)*EXV(-4.922_r8/T9) *RHO

                     !-B11 (ALPHA,N) N14                                 CF88
                     V(12)= 6.97D+12*TM23*EXV(-28.234_r8*TM13-(T9/0.140_r8)**2)*(1.0D0+0.015_r8 &
                           *TP13+8.115_r8*TP23+0.838_r8*T9+39.804_r8*TP43+10.456_r8*TP53)+1.79D+00 &
                           *TM32*EXV(-2.827_r8/T9)+1.71D+03*TM32*EXV(-5.178_r8/T9)+4.49D+06 &
                           *TP35*EXV(-8.596_r8/T9) *RHO

                     !-B11 (P,ALPHA) 1HE4                                NACRE
                     ! *** better to not use nacre fit, but tables
                     ! *** included by Falk. To switch to table.
                     V(13)=2.68D+12*TM23*EXV(-12.097_r8*TM13) *(1.D0 +1.62_r8*T9  &
                           -1.31_r8*T9*T9 +0.260_r8*T9*T9*T9) &
                           +2.12D+06*TM32*EXV(-1.724_r8/T9)*RHO*PROTOFAK
                  ENDIF
                  !------------------------------------------------------------
               END IF

               !------------------------------------------------------------
               !               CNO TRI CYCLE

               !------------------------------------------------------------
               IF(T9 .GT. 5.D-3) THEN
                  !**     IF(T9 .GT. 1.5D-2) THEN
                  !-C12(P,G)N13                          !CF88
                  V(14)=(2.04D+07/TP23*EXV(-13.690_r8/TP13-(T9/1.5_r8)**2)*(ONE+0.030_r8*TP13 &
                        +1.19_r8*TP23+0.254_r8*T9+2.06_r8*TP43+1.12_r8*TP53) + 1.08D+05/TP32 &
                        *EXV(-4.925_r8/T9) +2.15D+05/TP32*EXV(-18.179_r8/T9)) *RHO*PROTOFAK

                  ! *** NACRE
                  V(14)=(2.00D+07*TM23*exv(-13.692_r8*TM13- (T9/0.46_r8)**2)*(1.D0 &
                        +9.89_r8*T9-59.8D0*T9*T9 + 266.D0*T9*T9*T9)+1.00E+05*TM32 &
                        *exv(-4.913_r8/T9) + 4.24D+05*TM32* exv(-21.62_r8/T9))
                  V(14)=V(14)*RHO*PROTOFAK



                  !-N13(G,P)C12	REVERSE FORMULA FROM CF88
                  V(73) = V(14)/(RHO*PROTOFAK) * 8.84d+09*TP32 * exv(-22.553D0/T9)

                  !-C13(P,G)N14   CF1988
                  V(15)=(8.01D+07/TP23*EXV(-13.717_r8/TP13-(T9/2.0_r8)**2)*(1._r8+0.030_r8*TP13 &
                        +0.958_r8*TP23+0.204_r8*T9+1.39_r8*TP43+0.753_r8*TP53) &
                        + 1.21D+06/TP65  *EXV(-5.701_r8/T9))*RHO*PROTOFAK

                  ! *** NACRE
                  V(15)= (9.57D+07*TM23*exv(-13.720_r8*TM13-T9*T9) *(1+3.56_r8*T9) + &
                       1.50D+06*TM32*exv(-5.930_r8/T9) + exp(-12.057_r8/T9)* 6.83D+05/T9**0.864_r8)*&
                       RHO*PROTOFAK



                  ! ***  C-N13(P,G)O14..N14   CF88
                  V(16) = 4.04d+07*TM23*exv(-15.202_r8*TM13-(T9/1.191_r8)**2) &
                  *(1.0d0+0.027d0*TP13-0.803d0*TP23-0.154d0*T9+5.00d0*TP43 &
                  +2.44_r8*TP53)+2.43d+05*TM32*exv(-6.348d0/T9)*RHO*PROTOFAK

                  !-N14(P,G)O15..N15               CF 1988
                  V(17)=(4.90D7/TP23*EXV(-15.228_r8/TP13-(T9/3.294_r8)**2) &
                   * (ONE+0.027_r8*TP13-0.778_r8*TP23-0.149_r8*T9+0.261_r8*TP43+0.127_r8*TP53) &
                   + 2.37D3/TP32*EXV(-3.011_r8/T9)+2.19D4*EXV(-12.530_r8/T9)) &
                   * RHO * PROTOFAK


                  !    REACTION N14PG  IMBRIANI ET AL. 2005 (LUNA) (T9 < 2)
                  imb=1  !recommended
                  !        imb=2  !lower
                  !        imb=3  !upper
                  select case(imb)
                  case(1)
                     i1= 3.12d+7
                     i2= -15.193d0
                     i3= 0.486d0
                     i4= 0.782d0
                     i5= -1.50d0
                     i6= 17.97d0
                     i7= -3.32d0
                     i8= 2.11d+3
                     i9= -2.998d0
                     i10= 8.42d+2
                     i11= 0.0682d0
                     i12= -4.891d0
                  case(2)
                     i1= 2.76d+7
                     i2= -15.193d0
                     i3= 0.503d0
                     i4= 0.804d0
                     i5= -1.40d0
                     i6= 15.82d0
                     i7= -3.32d0
                     i8= 2.03d+3
                     i9= -2.998d0
                     i10= 8.44d+2
                     i11= 0.0682d0
                     i12= -4.987d0
                  case(3)
                     i1= 3.44d+7
                     i2= -15.193d0
                     i3= 0.475d0
                     i4= 0.771d0
                     i5= -1.59d0
                     i6= 19.83d0
                     i7= -3.30d0
                     i8= 2.18d+3
                     i9= -2.997d0
                     i10= 8.42d+2
                     i11= 0.0681d0
                     i12= -4.807d0
                  end select
                  if(t9.le.2._r8)then
                     v(17) = (i1/T9**TWO_THIRDS * exv(i2/T9**THIRD - (T9/i3)**2._r8) * &
                        (i4 + i5*T9 + i6*T9**2._r8 + i7*T9**3._r8) + i8/T9**(3._r8/2._r8) * exv(i9/T9) + &
                        i10*T9**i11 * exv(i12/T9))* RHO * protofak

                     !        call N14TEST(t9, rateout)
                     !        v(17) = rateout * RHO * PROTOFAK

                     !        else
                     !         print *,'WARNING: N14(p,g)_LUNA in VITAL is out of range'
                     !        print *,'rate used: CF88'
                  end if


                  !-N15(P,A)C12         FACT F1   FCZ1983 = CF88
                  F1=0.1_r8
                  V(18)=(1.08D+12/TP23*EXV(-15.251_r8/TP13-(T9/0.522_r8)**2)*(1._r8+0.027_r8 &
                   *TP13+2.62_r8*TP23+0.501_r8*T9+5.36_r8*TP43+2.6_r8*TP53) + 1.19D+08/TP32 &
                   *EXV(-3.676_r8/T9)+5.41D+08/TP12*EXV(-8.926_r8/T9) + F1 *(4.72D+08 &
                       /TP32*EXV(-7.721_r8/T9) +2.2D+09/TP32*EXV(-11.418_r8/T9)))* &
                    RHO*PROTOFAK
                  ! *** NACRE
                  V(18)=(1.12D12/T923*EXV(-15.253D0/T913-(T9/0.28D0)**2) &
                       *(1D0+4.95D0*T9+143D0*T9**2) &
                       +1.01D8/T932*EXV(-3.643D0/T9) &
                       +1.19D9/T932*EXV(-7.406D0/T9))* RHO * PROTOFAK

                  !-N15(P,G)O16      CF88
                  !      V(19)=(9.78D+08/TP23*EXV(-15.251_r8/TP13-(T9/0.45_r8)**2)*(ONE+0.027_r8
                  !     1*TP13+0.219_r8*TP23+0.042_r8*T9+6.83_r8*TP43+3.32_r8*TP53) + 1.11D+04/TP32
                  !     2*EXV(-3.328_r8/T9)+1.49D+04/TP32*EXV(-4.665_r8/T9)+3.8D+06/TP32
                  !     3*EXV(-11.048_r8/T9)) * RHO*PROTOFAK
                  ! *** NACRE
                  V(19)=(1.08D9/T923*EXV(-15.254D0/T913-(T9/0.34D0)**2) &
                       *(1D0+6.15D0*T9+16.4D0*T9**2) +9.23D3/T932*EXV(-3.597D0/T9) &
                       +3.27D6/T932*EXV(-11.024D0/T9))* RHO * PROTOFAK


                  !-O16(P,G)F17..O17    CF88
                  !      V(20)=(1.5D+08/(TP23*(ONE+2.13_r8*(ONE-EXV(-0.728_r8*TP23)))) *
                  !     1 EXV(-16.692_r8/TP13)) * RHO*PROTOFAK
                  ! *** NACRE
                  V(20)=(7.37D7/T9**0.82D0*EXV(-16.696D0/T913))*RHO*PROTOFAK

                  !-F17(G,P)O16..F17 --> O17	REVERSE FORMULA FROM CF88
                  V(74) = V(20)*3.03D+09*TP32*exp(-6.968D0/T9)/(RHO*PROTOFAK)

                  !-O17(P,A)N14        CF 88
                  !     F2=0.1_r8
                  !     F3=0.1_r8
                  !     V(21)=RHO*(1.53D+07/TP23*EXV(-16.712_r8/TP13-(T9/0.565_r8)**2)*(ONE+
                  !    1 0.025_r8*TP13+5.39_r8*TP23+0.940_r8*T9+13.5_r8*TP43+5.98_r8*TP53)+2.92D+06*T9*
                  !    2 EXV(-4.247_r8/T9) + F2 *(4.81D+10*T9*EXV(-16.712_r8/TP13-(T9/0.04_r8)**2)
                  !    3 +5.05D-05/TP32*EXV(-0.723_r8/T9))+F3 *1.31D+01/TP32*EXV(-1.961_r8/T9))
                  !    4 

                  !  Landre et al 1990 (A&A 240, 85)
                  F1=0.2_r8
                  F2=0.2_r8
                  V(21)=(1.53D+07/TP23*EXV(-16.712_r8/TP13-(T9/0.565_r8)**2)*(1._r8+ &
                   0.025_r8*TP13+5.39_r8*TP23+0.940_r8*T9+13.5_r8*TP43+5.98_r8*TP53)+2.92D+06*T9* &
                   EXV(-4.247_r8/T9) + 1.78D+05/TP23*EXV(-16.67_r8/TP13) &
                   /(0.479_r8*TP23+0.00312_r8)**2 + F1*2.8D+11*T9*EXV(-16.67_r8/TP13 &
                   -(T9/0.04_r8)**2) + F1 *2.94D-3/TP32*EXV(-0.767_r8/T9) &
                   + F2 * 98._r8/TP32*EXV(-2.077_r8/T9)) * RHO*PROTOFAK

                  !-O17(P,G)F18--> O18        FACT F3       FCZ1977   CF88
                  !      F3=F3
                  !     V(22)=(7.97D+07*T56/TP32*EXV(-16.712_r8/T13)+1.51D+08/TP23*EXV(-16.7_r8
                  !    112/TP13)*(ONE+0.025_r8*TP13-0.051_r8*TP23-8.82D-03*T9)+1.56D+05/T9*EXV
                  !    2( -6.272_r8/T9)+ F3 *1.31D+01/TP32*EXV(-1.961_r8/T9)) * RHO

                  !  Landre et al 1990 (A&A 240, 85)
                  F1=0.2_r8
                  F2=0.2_r8
                  TT=T9/(ONE+2.69_r8*T9)
                  T56=TT**FIVE_SIXTHS
                  T13=TT**THIRD
                  V(22)=(7.97D+07*T56/TP32*EXV(-16.712_r8/T13) &
                        + 1.51D+08/TP23 *EXV (-16.712_r8/TP13) &
                        * (1._r8+0.025_r8*TP13-0.051_r8*TP23-8.82D-3*T9) &
                        + 1.56D+05/T9 * EXV(-6.272_r8/T9) &
                        + F1 * 3.16D-5/TP32 *EXV(-0.767_r8/T9) &
                        + F2 * 98._r8/TP32 * EXV(-2.077_r8/T9)) * RHO*PROTOFAK

                  !-O18(P,A)N15         FACT F4    FCZ1985
                  !      F4=0.1_r8
                  !---------( F4=0. as in CF 88)
                  !       F4=0.
                  !      V(23)=(3.63D+11/TP23*EXV(-16.729_r8/TP13-(T9/1.361_r8)**2)*(ONE+0.025_r8
                  !     1 *TP13+1.88_r8*TP23+0.327_r8*T9+4.66_r8*TP43+2.06_r8*TP53)+2.35D-14/TP32*
                  !     2 EXV(-0.253_r8/T9)+2.66D+04/TP32*EXV(-1.67_r8/T9)+2.41D+09/TP32*EXV
                  !     3 (-7.638_r8/T9)+1.46D+09/T9*EXV(-8.31_r8/T9)+F4*(7.88D14/(TP23*(
                  !     4  0.439_r8*(ONE+5.18_r8*TP23)**2+0.561_r8))*EXV(-16.729_r8/TP13-0.534_r8*TP23)))
                  !     5    *RHO*PROTOFAK

                  ! *** O18(P,A)N15  CF88
                  V(23) = 3.63d+11/TP23*exv(-16.729_r8/TP13-(T9/1.361_r8)**2)*(1.0_r8+0.025_r8 &
                        *TP13+1.88_r8*TP23+0.327_r8*T9+4.66_r8*TP43+2.06_r8*TP53)+9.90d-14 &
                        /TP32*exv(-0.231_r8/T9)+2.66d+04/TP32*exv(-1.670_r8/T9)+2.41d+09 &
                        /TP32*exv(-7.638_r8/T9)+1.46d+09/T9*exv(-8.310_r8/T9) &
                      *RHO*PROTOFAK

                  !--------------------------------
                  !      IF(INETW .EQ. 1 .AND. IBURN .EQ. 1) THEN
                  !      DO LL=1,17
                  !      V(LL)=V(LL)*4.
                  !      END DO
                  !      END IF
                  !---------------------------------

               END IF


               !--------------------------------------------------
               !           Ne-- Na    and Mg-- Al cycles
               !--------------------------------------------------
               IF(INETW .EQ. 1 .AND. ICYCL .EQ. 1 .AND. T9 .GT. 5.D-3) THEN
                  !-O18(P,G)F19            CF1987-88
                  V(24)=(3.45E+08/TP23*EXV(-16.729_r8/TP13-(T9/0.139_r8)**2)*(1._r8+0.025_r8 &
                   *TP13+2.26_r8*TP23+0.394_r8*T9+30.56_r8*TP43+13.55_r8*TP53)+1.25E-15/TP32 &
                   *EXV(-0.231_r8/T9)+1.64E+02/TP32*EXV(-1.67_r8/T9)+1.28E+04*TP12 &
                   *EXV(-5.098_r8/T9)) * RHO *PROTOFAK
                  !------------ F ----------------------------------
                  !-F19(P,A)O16            CF1987-88
                  G1=ONE+FOUR*EXV(-2.09_r8/T9)+7._r8*EXV(-16.44_r8/T9)
                  V(25)=(3.55E+11/TP23*EXV(-18.113_r8/TP13-(T9/0.845_r8)**2)*(1._r8+0.023_r8 &
                   *TP13+1.96_r8*TP23+0.316_r8*T9+2.86_r8*TP43+1.17_r8*TP53)+3.67E+06/TP32 &
                   *EXV(-3.752_r8/T9)+3.07E+08*EXV(-6.019_r8/T9)) / G1 &
                   * RHO *PROTOFAK
                  !-F19(P,G)NE20           CF1987-88
                  G2=G1
                  V(26)=(6.04E+07/TP23*EXV(-18.113_r8/TP13-(T9/0.416_r8)**2)*(1._r8+0.023_r8 &
                   *TP13+2.06_r8*TP23+0.332_r8*T9+3.16_r8*TP43+1.30_r8*TP53)+6.32E+02_r8/TP32 &
                   *EXV(-3.752_r8/T9)+7.56E+04_r8/TP27*EXV(-5.722_r8/T9))/G2 &
                   * RHO*PROTOFAK

                  !------------- Ne 20  21  22-----------------------------
                  !=================== Art Champagne 1994 =========================
                  ! Ne 20 (p,g) Na 21     Champagne 1994
                  FA=T9**2*(ONE+0.0116_r8/TP23)**2
                  FC=-DSQRT(T9/0.21_r8)
                  FC=ONE+2.67_r8*EXP(FC)
                  V(27)=(2.05D8/TP23*EXV(-19.452_r8/TP13)*FC &
                        +9.55D6*EXV(-19.452_r8/TP13)/FA &
                        +18.0_r8/TP32*EXV(-4.257_r8/T9)+9.83_r8/TP32*EXV(-4.619_r8/T9) &
                        +1.59D5*T9**(1.49_r8)*EXV(-12.908_r8/T9)) &
                        * RHO *PROTOFAK
                  !
                  ! Ne 21 (p,g) Na 22     Champagne 1994
                  FNE21=0.10_r8
                  V(28)=(3.4D8/TP23*EXV(-19.41_r8/TP13)&
                             *(1._r8+0.56_r8*EXV(-(16.7_r8*T9-1._r8)**2))&
                             +6.12_r8/TP32*EXV(-1.403_r8/T9)+1.35D4/TP32*EXV(-3.008_r8/T9)&
                             +3.12D6*T9**(-0.72_r8)*EXV(-8.268_r8/T9**(+0.67_r8))&
                             + FNE21*1.1D-3/TP32*EXV(-1.114_r8/T9))&
                            * RHO * PROTOFAK
                  ! Ne 22 (p,g) Na 23     Champagne 1994
                  FNE22A=0.10_r8
                  FNE22B=0.10_r8
                  V(29)=(1.05D9/TP23*EXV(-19.43_r8/TP13)&
                             +1.24D-9/TP32*EXV(-0.414_r8/T9)+0.029_r8/TP32*EXV(-1.752_r8/T9)&
                             +9.30D4*T9**(-1.174_r8)*EXV(-5.100_r8/T9)&
                             +5.71D5*T9**(-0.249_r8)*EXV(-7.117_r8/T9)&
                             +FNE22A*3.25D-4/TP32*EXV(-0.789_r8/T9)&
                             +FNE22B*0.10_r8/TP32*EXV(-1.161_r8/T9))&
                             * RHO * PROTOFAK
                  !=================== Art Champagne 1994 =========================
                  !-NE20(P,G)NA21...NE21   CF1987-88
                  !      V(27)=(9.55E+06*EXV(-19.447_r8/TP13)/(T9**2*(ONE+0.0127_r8/TP23)**2)
                  !     & +2.05E+08/TP23*EXV(-19.447_r8/TP13)*(ONE+2.67_r8*EXV(-SQRT(T9/0.210_r8)))
                  !     & +1.80E+01/TP32*EXV(-4.242_r8/T9)+1.02E+01/TP32*EXV(-4.607_r8/T9)
                  !     & +3.60E+04/TP14*EXV(-11.249_r8/T9)) *RHO

                  !-NE21(P,G)NA22     CF1987-88
                  !      F4=0.1_r8               !  O TO 1 FACTOR
                  !      V(28)=(4.37E+08/TP23*EXV(-19.462_r8/TP13)+5.85_r8/TP32*EXV(-1.399_r8/T9)
                  !     & +1.29E+04/TP32*EXV(-3.009_r8/T9)+3.15E+05/TP35*EXV(-5.763_r8/T9)
                  !     & + F4 *(2.95E+08/TP23*EXV(-19.462_r8/TP13-(T9/0.058_r8)**2)*(ONE+0.021_r8
                  !     & *TP13+13.29_r8*TP23+1.99_r8*T9+124.1_r8*TP43+47.29_r8*TP53)+0.780_r8/TP32
                  !     & *EXV(-1.085_r8/T9)))* RHO
                  !-NE22(P,G)NA23            CF1987-88
                  !      F7=0.1_r8
                  !      V(29)=(1.15E+09/TP23*EXV(-19.475_r8/TP13) + 9.77E-12/TP32*EXV
                  !     &(-0.348_r8/T9) + 8.96E+03/TP32*EXV(-4.84_r8/T9) + 6.52E+04/TP32
                  !     &*EXV(-5.319_r8/T9) +   7.97E+05/TP12*EXV(-7.418_r8/T9) +  F7
                  !     &* 1.63E-01/TP32*EXV(-1.775_r8/T9)) * RHO
                  !-------------- NA  22  23  -----------------------------------
                  !C-NA22(B+)NE22       halflife=2.6088_r8 yr
                  !       V(25)=8.4188D-09
                  !-NA22(P,G)NA23            CF88
                  V(30)=9.63E-05_r8*TP32*EXV(-0.517_r8/T9)+2.51E+04_r8*T9*EXV(-2.013_r8/T9)
                  V(30)=V(30)*RHO*PROTOFAK
                  !=================== Art Champagne 1994 =========================
                  ! Na 23 (p,a) Na 20     Champagne 1994
                  FNA23=0.10_r8
                  V(31)=(1.26D10/TP23*EXV((-20.758_r8/TP13)-(T9/0.13_r8)**2)&
                             *(1._r8+0.02_r8*TP13-13.8_r8*TP23-1.93_r8*T9+234._r8*TP43+83.6_r8*TP53)&
                             +4.38_r8/TP32*EXV(-1.979_r8/T9)&
                             +6.50D6*T9**(-1.336_r8)*EXV(-6.490_r8/T9)&
                             +1.19D8*T9**(-1.055_r8)*EXV(-11.411_r8/T9)&
                             +FNA23*9.91D-14/TP32*EXV(-0.418_r8/T9))&
                             * RHO * PROTOFAK
                  ! Na 23 (p,g) Mg 24     Champagne 1994
                  FNA23P=0.10_r8
                  V(32)=(2.47D9/TP23*EXV(-20.758_r8/TP13)+91.9_r8/TP32*EXV(-2.789_r8/T9)&
                         +1.72D4/TP32*EXV(-3.433_r8/T9)&
                         +3.44D4*T9**(0.323_r8)*EXV(-5.219_r8/T9)&
                         +FNA23P*2.34D-4/TP32*EXV(-1.590_r8/T9))&
                         * RHO * PROTOFAK
                  !=================== Art Champagne 1994 =========================
                  !-NA23(P,A)NE20            CF1987-88                             FCZ1983
                  !      F8=0.1_r8               ! O TO 1 FACTOR
                  !      V(31)=(8.56E+09/TP23*EXV(-20.766_r8/TP13-(T9/0.131_r8)**2)*(ONE+0.02_r8*TP13
                  !     1 +8.21_r8*TP23+1.15_r8*T9+44.36_r8*TP43+15.84_r8*TP53)+4.02_r8/TP32*EXV(-1.99_r8/T9)
                  !     2 +1.18E+04/TP54*EXV(-3.148_r8/T9)+8.59E+05*TP43*EXV(-4.375_r8/T9) + F8
                  !     3 *3.06E-12/TP32*EXV(-0.447_r8/T9)) * RHO
                  !-NA23(P,G)MG24           CF1987-88
                  !      G3=ONE+1.5_r8*EXV(-5.105_r8/T9)
                  !      V(32)=(2.93E+08/TP23*EXV(-20.766_r8/TP13-(T9/0.297_r8)**2)*(ONE+0.02_r8*TP13
                  !     1 +1.61_r8*TP23+0.226_r8*T9+4.94_r8*TP43+1.76_r8*TP53)+93.4_r8/TP32*EXV(-2.789_r8/T9)
                  !     2+1.89E+04/TP32*EXV(-3.434_r8/T9)+5.10E+04*TP15*EXV(-5.51_r8/T9))/G3*RHO
                  !     & 
                  !-------------------- MG --------------------------------------
                  !-MG24(P,G)AL25...MG25    CF88 (see Trautvetter and Rolfs)
                  G4=ONE+5._r8*EXV(-15.882_r8/T9)
                  V(33)=(5.60E+08_r8/TP23*EXV(-22.019_r8/TP13)*(1._r8+0.019_r8*TP13-0.173_r8*TP23&
                   -0.023_r8*T9)+1480._r8/TP32*EXV(-2.484_r8/T9)+4000._r8*EXV(-4.18_r8/T9))/G4&
                     * RHO*PROTOFAK
                  !-MG25(P,G)AL26G (GROUND STATE)
                  !------( Art Champagne 1992) ----------
                  V(34)=(3.89D-8/TP32*EXV(-0.667_r8/T9)&
                        +2.06D-5/TP32*EXV(-1.070_r8/T9)&
                        +8.07D-2/TP32*EXV(-2.199_r8/T9)&
                        +4620._r8/TP32*EXV(-3.528_r8/T9)&
                        +2.17D+04*T9**(-0.075_r8)*EXV(-3.980_r8/T9))&
                        * RHO*PROTOFAK

                  !     VMG25PGAL26=3.57E+09/TP23*EXV(-22.031_r8/TP13-(T9/0.06_r8)**2.) *
                  !    1(ONE+0.019_r8*TP13+7.669_r8*TP23+1.015_r8*T9+167.4_r8*TP43+56.35_r8*TP53)
                  !    2+3.07E-13/TP32*EXV(-0.435_r8/T9)+1.94E-07/TP32*EXV(-0.673_r8/T9)
                  !    3+3.15E-05/T9**(3.4_r8)*EXV(-1.342_r8/T9-(T9/13.)**2)
                  !    4+1.77E+04*TP58*EXV(-3.049_r8/T9-(T9/13.)**2.)
                  !     V(26)= 0.8_r8 * VMG25PGAL26 * RHO     ! CF88
                  !----------------------------------------------------------------------
                  !     call SPLINT( TGRID , VMG25PG , Y2GMG25 , ngrid , T9 , VH )
                  !      Iliadis 90
                  !      V(29)=RHO * ( 10. ** VH )
                  !----------------------------------------------------------------------
                  !-MG25(P,G)AL26* (ISOMERIC STATE)  Mg26   ILiadis et al 90
                  !------( Art Champagne 1994) ----------
                  V(35)=(9.12D-9/TP32*EXV(-0.667_r8/T9)&
                        +4.27D-6/TP32*EXV(-1.070_r8/T9)&
                        +2.84D-2/TP32*EXV(-2.199_r8/T9)&
                        +690._r8/TP32*EXV(-3.528_r8/T9)&
                        +1.65D+04*T9**(-0.145_r8)*EXV(-4.436_r8/T9))&
                        * RHO*PROTOFAK

                  !     V(31)= 0.2_r8 * VMG25PGAL26 * RHO  ! CF88
                  !     call SPLINT( TGRID , VMG25PI , Y2IMG25 , ngrid , T9 , VH )
                  !     V(30)= RHO * ( 10. ** VH )

                  !-MG26(P,G)AL27
                  !---- Champagne etal + Iiadis etal + Buchmann et al + Endt
                  F1=0.10_r8
                  V(36)=(8.54D-12/TP32*EXV(-0.605_r8/T9)&
                        +2.747D-6/TP32*EXV(-1.219_r8/T9)&
                        +0.0129_r8/TP32*EXV(-1.728_r8/T9)&
                        +8.06_r8/TP32*EXV(-2.537_r8/T9)&
                        +1450._r8/TP32*EXV(-3.266_r8/T9)&
                        +4.03D+04/TP32*EXV(-3.784_r8/T9)&
                        +8.82D+04*T9**(-0.21_r8)*EXV(-4.194_r8/T9)&
                        + F1 * 1.93D-5/TP32 * EXV(-1.044_r8/T9))&
                        * RHO*PROTOFAK

                  !     call SPLINT( TGRID , VMG26PG , Y2MG26 , ngrid , T9 , VL )
                  !      V(32)=RHO * ( 10. ** VL )
                  !-------------------- AL ------------------------------------
                  !-AL26(P,G)SI27...AL27            VOGELAAR89
                  IF(T9 .LT. 0.1_r8) THEN
                     !.... lower limit according to Art (1994) ....
                     V(37)=(4.67D-10/TP32*EXV(-0.789_r8/T9)+3.70D-8/TP32*EXV(-1.029_r8/T9))&
                           * RHO *PROTOFAK
                  ELSE
                     V(37)=(8.97_r8/TP32 *EXV(-2.191_r8/T9) + 473._r8/TP32*EXV(-3.22_r8/T9)&
                            +7763._r8/T9*EXV(-3.944_r8/T9)) * RHO*PROTOFAK
                  END IF
                  !............... old one .below T9 .lt. 0.1..................
                  !      V(37)=( 8.744E+14/TP23*EXV(-27.602_r8/TP13*(ONE+.158*T9-.0107*T9*T9+
                  !     &     3.716E-4*T9*T9*T9)) + 3.26E-10/TP32*EXV(-.805/T9) + 3.26E-3
                  !     &     /TP32*EXV(-1.453_r8/T9) + 8.97_r8/TP32*EXV(-2.191_r8/T9) + 473./TP32
                  !     &     *EXV(-3.22_r8/T9) + 7763./T9*EXV(-3.944_r8/T9) )* RHO

                  !C-AL26G(B+,G)MG26
                  !      V(34)=3.06E-14
                  !C-AL26*(B+,G)MG26
                  !***      V(34)=1.6E-01
                  !C
                  !-AL27(P,A)MG24              CHAMPAGNE 1994
                  ARG1=-23.25_r8/TP13 - 3.57_r8*T9**2
                  ARG2=-2.148_r8/T9**(1.293_r8)
                  F27A=0.17_r8
                  F27B=0.10_r8
                  V(38)=(4.71D+05/TP23*EXV(ARG1) * (1._r8+0.018_r8*TP13-7.29_r8*TP23&
                        -0.914_r8*T9+77.20_r8*TP43+24.60_r8*TP53)&
                        +2.23D+4*T9**(3.989_r8)*EXV(ARG2)&
                        +F27A*1.29D-9/TP32*EXV(-0.836_r8/T9)&
                        +F27B*2.73D-3/TP32*EXV(-2.269_r8/T9)&
                        +F27B*0.026_r8/TP32*EXV(-2.492_r8/T9))&
                        * RHO * PROTOFAK
                  !
                  !      SIGMAV = ( 1.82E-10*EXV(-.853/T9) + 6.6E-11*EXV(-ONE/T9)
                  !     &       + 2.81E-3*EXV(-2.27_r8/T9) + 2.64E-2*EXV(-2.51_r8/T9) + .198
                  !     &       *EXV(-3.3_r8/T9) + .528*EXV(-3.67_r8/T9) + 2.81_r8*EXV(-4.55_r8/T9)
                  !     &       + 495.*EXV(-7.34_r8/T9) ) / TP32
                  !      V(38)=RHO * SIGMAV
                  !
                  !-AL27(P,G)SI28             Champagne 1994
                  F17=0.17_r8
                  V(39)=(1.32D+9/TP23*EXV(-23.26_r8/TP13)&
                         +F17*3.22D-10/TP32*EXV(-0.836_r8/T9)&
                         +1.74_r8/TP32*EXV(-2.269_r8/T9) + 9.92_r8/TP32*EXV(-2.492_r8/T9)&
                         +42.9_r8/TP32*EXV(-3.273_r8/T9) + 134._r8/TP32*EXV(-3.654_r8/T9)&
                         +1.77D+4*T9**(0.53_r8)*EXV(-4.588_r8/T9))&
                         * RHO * PROTOFAK
                  !
                  !-AL27(P,G)SI28                   CF 1988
                  !      F13=0.1_r8
                  !      G7=ONE+EXV(-9.792_r8/T9)/3.+2.*EXV(-11.773_r8/T9)/3.
                  !      V(39)=(1.67E+08/TP23*EXV(-23.261_r8/TP13-(T9/0.155_r8)**2)*(ONE+0.018_r8
                  !     1 *TP13+5.81_r8*TP23+0.728_r8*T9+27.31_r8*TP43+8.71_r8*TP53)+2.2_r8/TP32
                  !     2 *EXV(-2.269_r8/T9)+12.2_r8/TP32*EXV(-2.491_r8/T9)+1.5E+04*T9*EXV(-4.112_r8
                  !     3 /T9) + F13*(6.5E-10/TP32*EXV(-0.853_r8/T9)+1.63E-10/TP32
                  !     4 *EXV(-1.001_r8/T9)))/ G7 * RHO 

                  !-Si28(p,G)P29 (B+)Si29         CF88
                  V(40)=(1.64E+8/TP23*EXV(-24.449_r8/TP13-(T9/2.91_r8)**2)*(ONE+0.017_r8&
                        *TP13-4.11_r8*TP23-0.491_r8*T9+5.22_r8*TP43+1.58_r8*TP53)&
                        +3.52E+02/TP32*EXV(-4.152_r8/T9)+6.30E+05/TP32*&
                        EXV(-18.505_r8/T9)+1.69E+03*EXV(-14.518_r8/T9))&
                        *RHO*PROTOFAK
                  !-Si29(p,G)P30(B+)SI30           CF88
                  V(41)=(3.26E+09/TP23*EXV(-24.459_r8/TP13-(T9/0.256_r8)**2)&
                         *(ONE+0.017_r8*TP13+4.27_r8*TP23+0.509_r8*T9+15.40_r8*TP43+4.76_r8*&
                         TP53)+2.98E+03/TP32*EXV(-3.667_r8/T9)+3.94E+04/TP32*&
                         EXV(-4.665_r8/T9)+2.08E+04*TP12*EXV(-8.657_r8/T9))&
                         *RHO*PROTOFAK
                  !-Si30(p,G)P31                   CF88
                  V(42)=(4.25E+08/TP23*EXV(-24.468_r8/TP13-(T9/0.670_r8)**2)&
                         *(ONE+0.017_r8*TP13+0.15_r8*TP23+0.018_r8*T9+5.53_r8*TP43+1.68_r8*&
                         TP53)+1.86E+04/TP32*EXV(-5.601_r8/T9)+3.15E+05/TP32*&
                         EXV(-6.961_r8/T9)+2.75E+05/TP12*EXV(-10.0672_r8/T9))&
                         *RHO*PROTOFAK
               END IF

               !----------------------------------------------------------------
               !               HELIUM BURNING
               !----------------------------------------------------------------
               !--- HE 4
               !-HE4(2A,G)C12                                   FCZ1983 = CF88
               IF(T9.GT.0.08_r8)THEN
                  F14=0.1_r8
                  V(43)=(2.79D-08/T9**3*EXV(-4.4027_r8/T9)+F14*1.35D-07 *TM32* &
                    EXV(-24.811_r8/T9))*RHO**2
               ELSE IF(T9.GE.5.E-3) THEN
                  F14=0.1_r8
                  V1=7.4D+05/TP32*EXV(-1.0663_r8/T9)+4.164D+09/TP23*EXV(-13.490_r8/TP13&
                      -(T9/0.098_r8)**2)*(ONE+0.031_r8*TP13+8.009_r8*TP23+1.732_r8*T9+49.883_r8*&
                      TP43+27.426_r8*TP53)
                  V2=1.30D+02/TP32*EXV(-3.3364_r8/T9)+2.510D+07/TP23*EXV(-23.570_r8/TP13&
                      -(T9/0.235_r8)**2)*(ONE+0.018_r8*TP13+5.249_r8*TP23+0.650_r8*T9+19.176_r8*TP43&
                      +6.034_r8*TP53)
                  V(43)=(2.9D-16*V1*V2*(0.01_r8+0.2_r8*(ONE+4.*EXV(-(0.025_r8/T9)**3.263_r8))&
                   /(ONE+4.*EXV(-(T9/0.025_r8)**9.227_r8)))+F14*1.35D-07/TP32*EXV(-24.811_r8&
                     /T9))*(RHO**2)
               ELSE IF(T9.lt.5.E-3) THEN
                  V(43) = 0.d0
               END IF
               !        print *,'old V(43) * rho',V(43)
               !-HE4(2A,G)C12                                                NACRE REC
               ! *** better to not use nacre fit, but tables
               ! *** included by Falk. To switch to table.
               HE4HE4 = 2.43D+9 * TM23 * &
                   EXP(-13.490d0 * TM13 - (T9/0.15D0)**2) * &
                   (1.D0 + 74.5D0 * T9) + 6.09D+5 * TM32 * &
                   EXP(-1.054D0/T9)
               HE4BE8 = 2.76D+7 * TM23 * EXP(-23.570D0 * TM13 - (T9/0.4D0)**2.D0) *&
                   (1.D0 + 5.47D0 * T9 + 326.D0 * T9**2) +&
                   130.7D0 * TM32 * EXP(-3.338D0/T9) + 2.51D0 * TM32 * EXP(-20.307D0/T9)
               if (T9.le.0.03d0)then
                  V(43) = HE4HE4 * HE4BE8 * 3.07D-16 * (ONE - 29.1D0 * T9 + 1308.D0 * T9**2)
               else if (T9.gt.0.03d0)then
                  V(43) = HE4HE4 * HE4BE8 * 3.44D-16 * (ONE + 0.0158D0 * T9**(-0.65D0))
               end if


               !      HE4(2A,G)C12    reaclib JINA - Fynbo et al. 2005 Nature 433, 136-139
               if (t9 >= 0.01_r8) then
                  ! non resonant
                  a0r = -9.710520d-01
                  a1r =  0.000000d+00
                  a2r = -3.706000d+01
                  a3r =  2.934930d+01
                  a4r = -1.155070d+02
                  a5r = -1.000000d+01
                  a6r = -1.333330d+00

                  Vn77 = a0r + a1r/t9 + a2r/TP13 + a3r*TP13 + a4r*t9 + a5r*TP53 + a6r*lnt9
                  ! *** resonant 1
                  a0r = -2.435050d+01
                  a1r = -4.126560d+00
                  a2r = -1.349000d+01
                  a3r =  2.142590d+01
                  a4r = -1.347690d+00
                  a5r =  8.798160d-02
                  a6r = -1.316530d+01

                  Vr177 = a0r + a1r/t9 + a2r/TP13 + a3r*TP13 + a4r*t9 + a5r*TP53 + a6r*lnt9
                  ! *** resonant 2
                  a0r = -1.178840d+01
                  a1r = -1.024460d+00
                  a2r = -2.357000d+01
                  a3r =  2.048860d+01
                  a4r = -1.298820d+01
                  a5r = -2.000000d+01
                  a6r = -2.166670d+00

                  Vr277 = a0r + a1r/t9 + a2r/TP13 + a3r*TP13 + a4r*t9 + a5r*TP53 + a6r*lnt9

                  V(43) = (exp(Vn77) + exp(Vr177) + exp(Vr277)) * rho * rho
               else
                  V(43) = ZERO
               end if
               ! *** here I include the reverse rates of 3alpha: C12-->3alpha
               ! *** from public_torch. source of reverse rate: fxtsource

               v(80) = v(43)/(rho*rho)*2.00d+20*(t9**3)*exp(-84.424_r8/t9)


               !------- C 12
               IF(T9 .GE. 0.04D0) THEN
                  !-C12(A,G)O16                           CFHZ1985
!                       V(44)=(2.78D+08/T9**2./(ONE+0.0489_r8/TP23)**2.*EXV(-32.12_r8/TP13- &
!                      (T9/3.496_r8)**2.)+4.68D+08/T9**2./(ONE+0.2654_r8/TP23)**2.*EXV(-32.12_r8 &
!                      /TP13)+1.25D+03/TP32*EXV(-27.499_r8/T9)+1.43D-02*T9**5.*EXV(-15.541_r8 &
!                      /T9))*RHO
                  ! ===>C12(A,G) O16     CF 1988: <===
                  V(44)=(1.04e8_r8/T9**2/(ONE + 0.0489_r8/TP23)**2 * EXV(-32.12_r8/TP13 - (T9/3.496_r8)**2) &
                       +1.76e8_r8 / T9**2 / (ONE + 0.2654_r8/TP23)**2 * EXV(-32.120_r8/TP13) &
                       +1.25e3_r8/TP32 * EXV(-27.499_r8/T9) + 1.43e-2_r8*T9**5 * EXV(-15.541_r8/T9)) * RHO * 1.7_r8


                  !    REAZIONE C12AG  Buchmann1996  0.03_r8 < T9 < 2
                  ! *** corrected according to Buchmann1997.
                  Buch=1   !recommended rate
                  !         Buch=2   !upper limit
                  !         Buch=3   !lower limit
                  if(Buch == 1)then
                     p1= -3.5738d+7
                     p2= 0.65711d0
                     p3= 30.834d0
                     p4= 0.40104d0
                     p5= 5.0464d+8
                     p6= 0.66456d0
                     p7= 31.332d0
                     p8= -1.6054d0
                     p9= 16.272d0
                     p10=85028.d0
                     p11=8844.5d0
                     p12=33.033d0
                  else if(Buch == 2)then
                     p1= 1.7399e+8_r8
                     p2= 13.608d0
                     p3= 40.171d0
                     p4= -11.013d0
                     p5= 7.6099d+8
                     p6= 0.47121d0
                     p7= 31.666d0
                     p8= -0.78718d0
                     p9= 14781d0
                     p10=56546.d0
                     p11=5881.1d0
                     p12=31.332d0
                  else if(Buch == 3)then
                     p1= -5.1029d+8
                     p2= 16.686d0
                     p3= 29.877d0
                     p4= 86.063d0
                     p5= 3.4183d+8
                     p6= 1.4064d0
                     p7= 30.955d0
                     p8= -0.71494d0
                     p9= 15.514d0
                     p10=1.7916d+5
                     p11=18617.d0
                     p12=36.057d0
                  end if
                  !        if(t9.ge.0.03_r8.and.t9.le.2.)then
                  !         V(44) = (p1/(T9**2.d0 * (1.d0 + p2/T9**(2.d0/3.d0))**2.) *
                  !     &    dexp((-p3/T9**(1.d0/3.d0)) - (T9**2.d0/p4)) +
                  !     &    p5/(T9**2.d0*(1.0d0 + p6/T9**(2.d0/3.d0))**2.) *
                  !     &    dexp(-p7/T9**(1.d0/3.d0)) +
                  !     &    p8/T9**(3.d0/2.d0)*dexp(-p9/T9) +
                  !     &    p10/T9**(2.d0/3.d0)*(1.d0 + p11*T9**(1.d0/3.d0)) *
                  !     &    dexp(-p12/T9**(1.d0/3.d0)))
                  !     &    * RHO
                  !        else
                  !         stop'Buch96 used in VITAL for C12agO16. Out of T range'
                  !        end if

                  !    REAZIONE C12AG  kunz02   (0.001_r8 < T9 < 10)
                  kunz=1  !recommended
                  !        kunz=2  !upper
                  !        kunz=3  !lower
                  if (kunz == 1)then
                     a0= 1.21d+8
                     a1= 6.06d-2
                     a2= 32.12d0
                     a3= 1.7d0
                     a4= 7.4d+8
                     a5= 0.47d0
                     a6= 32.12d0
                     a7= ZERO
                     a8= ZERO
                     a9= 1.53d+4
                     a10=2.d+6
                     a11=38.534d0
                     !A9= 3.06d+10
                  else if(kunz == 2)then
                     a0= 1.35d+8
                     a1= 5.45d-2
                     a2= 32.12d0
                     a3= 10.d0
                     a4= 9.4d+8
                     a5= 0.41d0
                     a6= 32.12d0
                     a7= ZERO
                     a8= ZERO
                     a9= 1.7d+4
                     a10=2.22d+6
                     a11=38.6d0
                     !A9= 3.77d+10
                  else if(kunz == 3)then
                     a0= 3.2d+7
                     a1= 3.5d-2
                     a2= 32.12d0
                     a3= 8.d-2
                     a4= 4.6d+8
                     a5= 0.262d0
                     a6= 32.12d0
                     a7= ZERO
                     a8= ZERO
                     a9= 1.4d+4
                     a10=1.9d+6
                     a11=38.67d0
                     !A9= 2.66d+10
                  end if
                  v(44) = (a0/(T9**2*(ONE + a1/TP23)**2) * exp((-a2/TP13) - (T9/a3)**2) + &
                      a4/(T9**2 * (ONE + a5/TP23)**2) * exp(-a6/TP13) + a7/TP32*exp(-a8*T9inv) +&
                      a9/TP23*(ONE + a10*TP13) * exp(-a11/TP13)) * RHO


                  !---- modification for c12(a,g) ------J. DeBoer+2016-----!
                  !                                            by MP       !
                  if ( c12ag_jdb2016 ) then
                     deboer = 2 ! 2=recomm. rate; 3=upper limit; 1= lower limit
                     rt_c12ag(1:ndat_c12ag-1) = rt_c12ag_in(1:ndat_c12ag-1,deboer)
                     rt_c12ag(1:ndat_c12ag-1) = log10( rt_c12ag(1:ndat_c12ag-1) )
                     ! ^_^ should just be done with standard 1d interp, clipping t9 to bounds...
                     if      ( t9 <= t9_in(1) ) then
                        i_c12 = 1
                        rate_c12ag_nd = 10._r8 ** rt_c12ag(i_c12) 
                     else if ( t9 >= t9_in(ndat_c12ag-1) ) then
                        i_c12 = ndat_c12ag-1
                        rate_c12ag_nd = 10._r8 ** rt_c12ag(i_c12) 
                     else
                        do idt = 1, ndat_c12ag - 2
                           if ( t9_in(idt) <= t9 .and. t9 < t9_in(idt+1) ) then
                              i_c12 = idt
                              h_intpl = ( t9 - t9_in(i_c12) ) / ( t9_in(i_c12 + 1) - t9_in(i_c12) )
                              ! *** c12(a,g)
                              rate_c12ag_nd = rt_c12ag(i_c12) + h_intpl * &
                                    ( rt_c12ag(i_c12 + 1) - rt_c12ag(i_c12) ) 
                              rate_c12ag_nd = 10._r8 ** rate_c12ag_nd
                              exit
                           end if
                        end do
                     end if
                     v(44) = rate_c12ag_nd * RHO
                  end if

                  !-------  C 13
                  !- 1 C  13 (1 ALPHA, 1 NEUT)  1 O 16           CF 1988
                  V(45)=(6.77D+15/TP23*EXV(-32.329_r8/TP13-(T9/1.284_r8)**2)&
                        *(ONE+0.013_r8*TP13+2.04_r8*TP23+0.184_r8*T9)&
                        +3.82D+5/TP32*EXV(-9.373_r8/T9)+1.41D+6/TP32*EXV(-11.873_r8/T9)&
                        +2.0D+9/TP32*EXV(-20.409_r8/T9)+2.92D+9/TP32*EXV(-29.283_r8/T9))
                  V(45)=V(45) *RHO

                  !      if (T9.lt.4.) then          ! NACRE
                  !      V(45)=(3.78D14/(T9*T9)*exv(-32.333_r8*TM13-(T9/0.71_r8)**2)*(1.D0+
                  !     & 46.8D0*T9- 292.D0*T9*T9+738.*T9*T9*T9)+2.30D+07*T9**0.45_r8*
                  !     & exv(-13.03_r8/T9)) *RHO
                  !      else
                  !      stop 'F: C13(alpha,n)O16 rate from NACRE compilation
                  !     & out of range'
                  !      end if

                  ! *** c13(a,n) Heil et al. 2008 Phys Rev C + priv comm. Heil 2009
                  ! *** note that the rate is tuned between t9 = 0.04_r8 and t9 = 1.
                  ! *** it is not safe to use this formula out of this temperature range.
                  ! *** e.g. at T9=10 this formula is giving negative rates....
                  ! *** so, use CF88 rate outside the allowed temperature range
                  ! *** (at low temperatures the fit is just giving very low rates,
                  ! *** correctly)
                  if (t9.le.1.0_r8) then
                     V(45)=(4.5439D14/(T9*T9)*exv(-32.48235_r8*TM13-(T9/1.63585_r8)**2)&
                      * (1.D0 + 20.05918D0*T9 - 71.98345D0*T9*T9 + 258.58608_r8*T9*T9*T9)&
                      - 5.247020661485D+07*T9**0.45_r8*exv(-12.00061_r8/T9))
                     V(45)=V(45) *RHO
                  end if

               END IF
               !
               IF(INETW .EQ. 1.AND. T9 .GE. 0.08_r8) THEN
                  !-N14(A,G)F18                                                 FCZ1975
                  V(46)=(7.78D+09*TM23*EXV(-36.031_r8*TM13-(T9/0.881_r8)**2)&
                    *(ONE+0.012_r8*TP13+1.45_r8*TP23+0.117_r8*T9+1.97_r8*TP43+0.406_r8*TP53)&
                    +2.36D-10*TM32*EXV(-2.798_r8/T9)+2.03D+00*TM32*EXV(-5.054_r8/T9)&
                    +1.15D+04*TM23*EXV(-12.31_r8/T9))*RHO
               END IF
               IF(T9 .GE.0.08_r8) THEN
                  ! O16(A,G)NE20                                                 CF 1988
                  !      V(47)=(9.37D09/TP23*EXV(-39.757_r8/TP13-(T9/1.586_r8)**2)+6.21D+1/TP32
                  !     &     *EXV(-10.297_r8/T9)+5.38D02/TP32*EXV(-12.226_r8/T9)
                  !     &     +13.0_r8*T9**2*EXV(-20.093_r8/T9))*RHO
                  ! O16(A,G)NE20                                                 NACRE
                  r_o16ag =&
                       (2.68D10/TP23*EXV(-39.760_r8/TP13-(T9/1.6d0)**2)+5.11D+1/TP32&
                       *EXV(-10.32d0/T9)+6.161D02/TP32*EXV(-12.200_r8/T9)&
                       +0.41d0*T9**2.966_r8*EXV(-11.900d0/T9))
                  v(47) = r_o16ag * RHO
               END IF

               !
               IF(INETW .EQ. 1 .AND. T9 .GE. 0.08_r8) THEN
                  IF(T9 .LT. 0.80D0) THEN
                     !-O18(A,G)NE22         Giesen et al. (1994)
                     CALL CUBIC(T9,VO18)
                     V(48)= VO18 * RHO
                  ELSE
                     !-O18(A,G)NE22          CF 1988
                     V(48)=(1.82D12/TP23*EXV(-40.057_r8/TP13-(T9/0.343_r8)**2)&
                           *(ONE+0.010_r8*TP13+0.988_r8*TP23+0.072_r8*T9+3.17_r8*TP43&
                           +0.586_r8*TP53)+7.54D00/TP32*EXV(-6.228_r8/T9)+3.48D01&
                           /TP32*EXV(-7.301_r8/T9)+6.23D03*T9*EXV(-16.987_r8/T9)&
                           +0.1_r8*(1.D-11/TP32*EXV(-1.994_r8/T9)))*RHO
                  END IF
               END IF
               !
               IF(T9 .GE. 0.08_r8)THEN
                  !-NE20(A,G)MG24
                  V(49)=(4.11D11/TP23*EXV(-46.766_r8/TP13-(T9/2.219_r8)**2)&
                       *(ONE+0.009_r8*TP13+0.882_r8*TP23+0.055_r8*T9+0.749_r8*TP43+0.119_r8&
                       *TP53)+5.27D03/TP32*EXV(-15.869_r8/T9)+6.51D03*TP12*&
                       EXV(-16.223_r8/T9)+0.1_r8*(4.21D01/TP32*EXV(-9.115_r8/T9)&
                       +3.2D01/TP23&
                       *EXV(-9.383_r8/T9)))/(ONE+5.*EXV(-18.960_r8/T9))*RHO
               END IF

               IF(INETW .EQ. 1) THEN
                  !- 1 NE 21 ( 1 ALPHA, 1 NEUT )  1 MG 24                         CF  1987
                  T9A   = T9 / (ONE + 0.0537_r8 * T9)
                  TT13  = T9A ** THIRD
                  TT56  = T9A ** FIVE_SIXTHS
                  GT9   = ONE + 1.5_r8 * EXV(-4.068_r8/T9)+2._r8 * EXV(-20.258_r8/T9)
                  V(50) = (4.94D+19*TT56 /TP32*EXV(-46.89_r8/TT13) + 2.66D+07/TP32 &
                         *EXV(-22.049_r8/T9) )* RHO /GT9
               END IF


               !---- modification for Ne22(a,n)/(a,g) -----------------------!
               !                                                by NN+MP     !
               ! ^_^ TODO: THE FILE SHOULD NOT BE READ EVERY PHYSICS CALL!!!!!!!!
               if ( ne22_michael ) then
                  ne22an = 1 ! 1=recomm. rate; 2=upper limit; 3= upper upper limit
                  ne22ag = 1 ! 1=recomm. rate; 2=upper limit; 3= upper upper limit
                  open(123, file = '../NPDATA/rate_ne22.dat', action = 'read')
                  do idt = 1, ndat_ne22
                     read(123,*) t9_in(idt), rt_ne22an_in(idt,1:3), rt_ne22ag_in(idt,1:3)
                  end do
                  close(123)

                  rt_ne22an(1:ndat_ne22) = rt_ne22an_in(1:ndat_ne22,ne22an)
                  rt_ne22ag(1:ndat_ne22) = rt_ne22ag_in(1:ndat_ne22,ne22ag)

                  rt_ne22an(1:ndat_ne22) = log10( rt_ne22an(1:ndat_ne22) )
                  rt_ne22ag(1:ndat_ne22) = log10( rt_ne22ag(1:ndat_ne22) )

                  ! ^_^ TODO: again, this should just be done with a standard 1d table interp, clipping t9 to table bounds
                  if ( t9 <= t9_in(1) ) then
                     i_ne22 = 1
                     rate_ne22ag_mw = 10**rt_ne22ag(i_ne22) 
                     rate_ne22an_mw = 10**rt_ne22an(i_ne22) 
                  else if ( t9 >= t9_in(ndat_ne22) ) then
                     i_ne22 = ndat_ne22
                     rate_ne22ag_mw = 10**rt_ne22ag(i_ne22) 
                     rate_ne22an_mw = 10**rt_ne22an(i_ne22)
                  else
                     do idt = 1, ndat_ne22 - 1
                        if ( t9_in(idt) <= t9 .and. t9 < t9_in(idt+1) ) then
                           i_ne22 = idt
                           h_intpl = ( t9 - t9_in(i_ne22) ) / ( t9_in(i_ne22 + 1) - t9_in(i_ne22) )
                           ! *** Ne22(a,g)
                           rate_ne22ag_mw = rt_ne22ag(i_ne22) + h_intpl * ( rt_ne22ag(i_ne22 + 1) - rt_ne22ag(i_ne22) ) 
                           rate_ne22ag_mw = 10._r8 ** rate_ne22ag_mw
                           ! *** Ne22(a,n)
                           rate_ne22an_mw = rt_ne22an(i_ne22) + h_intpl * ( rt_ne22an(i_ne22 + 1) - rt_ne22an(i_ne22) ) 
                           rate_ne22an_mw = 10._r8 ** rate_ne22an_mw
                           exit
                        end if
                     end do
                  end if
               end if

               !                                                             !
               !---- modification for Ne22(a,n)/(a,g) -----------------------!


               !---- modification for Ne22(a,n)/(a,g) ---Longland+2012-------!
               !                                                by MP        !
               ! ^_^ TODO: THE FILE SHOULD NOT BE READ EVERY PHYSICS CALL!!!!!!!!
               if ( ne22_longland ) then
                  ne22an = 3 ! 2=recomm. rate; 1=lower limit; 3= upper limit
                  ne22ag = 1 ! 2=recomm. rate; 1=lower limit; 3= upper limit
                  open(124, file = '../NPDATA/ag_lo12.dat', action = 'read')
                  read(124,*)
                  read(124,*)
                  read(124,*)
                  do idt = 1, ndat_ne22
                     read(124,*) t9_in(idt), rt_ne22ag_in(idt,1:3)
                  end do
                  close(124)
                  open(125, file = '../NPDATA/an_lo12.dat', action = 'read')
                  read(125,*)
                  read(125,*)
                  read(125,*)
                  do idt = 1, ndat_ne22
                     read(125,*) t9_in(idt), rt_ne22an_in(idt,1:3)
                  end do
                  close(125)

                  rt_ne22an(1:ndat_ne22) = rt_ne22an_in(1:ndat_ne22,ne22an)
                  rt_ne22ag(1:ndat_ne22) = rt_ne22ag_in(1:ndat_ne22,ne22ag)

                  rt_ne22an(1:ndat_ne22) = log10( rt_ne22an(1:ndat_ne22) )
                  rt_ne22ag(1:ndat_ne22) = log10( rt_ne22ag(1:ndat_ne22) )

                  ! ^_^ TODO: again, this should just be done with a standard 1d table interp, clipping t9 to table bounds
                  if ( t9 <= t9_in(1) ) then
                     i_ne22 = 1
                     rate_ne22ag_lo = 10._r8 ** rt_ne22ag(i_ne22) 
                     rate_ne22an_lo = 10**rt_ne22an(i_ne22) 
                  else if ( t9 >= t9_in(ndat_ne22) ) then
                     i_ne22 = ndat_ne22
                     rate_ne22ag_lo = 10**rt_ne22ag(i_ne22)
                     rate_ne22an_lo = 10**rt_ne22an(i_ne22)
                  else
                     do idt = 1, ndat_ne22 - 1
                        if ( t9_in(idt) <= t9 .and. t9 < t9_in(idt+1) ) then
                           i_ne22 = idt
                           h_intpl = ( t9 - t9_in(i_ne22) ) & 
                                /( t9_in(i_ne22 + 1) - t9_in(i_ne22) )
                           ! *** Ne22(a,g)
                           rate_ne22ag_lo = rt_ne22ag(i_ne22)  &
                                   + h_intpl *( rt_ne22ag(i_ne22 + 1) - rt_ne22ag(i_ne22) ) 
                           rate_ne22ag_lo = 10.d0**rate_ne22ag_lo
                           ! *** Ne22(a,n)
                           rate_ne22an_lo = rt_ne22an(i_ne22)  &
                                   + h_intpl *( rt_ne22an(i_ne22 + 1) - rt_ne22an(i_ne22) ) 
                           rate_ne22an_lo = 10.d0**rate_ne22an_lo
                           exit
                        end if
                     end do
                  end if
               end if

               !                                                             !
               !---- modification for Ne22(a,n)/(a,g) -----------------------!

               !---- modification for Ne22(a,n)/(a,g) ------Talwar+2015-----!
               !                                                by MP       !
               ! ^_^ TODO: THE FILE SHOULD NOT BE READ EVERY PHYSICS CALL!!!!!!!!
               if ( ne22_nd2015 ) then
                  ne22an = 3 ! 2=recomm. rate; 3=upper limit; 1= lower limit
                  ne22ag = 1 ! 2=recomm. rate; 3=upper limit; 1= lower limit
                  open(125, file = '../NPDATA/ne22a_mw15.dat', action = 'read')
                  read(125,*)
                  read(125,*)
                  read(125,*)
                  read(125,*)
                  do idt = 1, ndat_ne22-1
                     read(125,*) t9_in(idt), rt_ne22ag_in(idt,1:3), rt_ne22an_in(idt,1:3)
                  end do
                  close(125)

                  rt_ne22an(1:ndat_ne22-1) = rt_ne22an_in(1:ndat_ne22-1,ne22an)
                  rt_ne22ag(1:ndat_ne22-1) = rt_ne22ag_in(1:ndat_ne22-1,ne22ag)

                  rt_ne22an(1:ndat_ne22-1) = log10( rt_ne22an(1:ndat_ne22-1) )
                  rt_ne22ag(1:ndat_ne22-1) = log10( rt_ne22ag(1:ndat_ne22-1) )

                  ! ^_^ TODO: again, this should just be done with a standard 1d table interp, clipping t9 to table bounds
                  if ( t9 <= t9_in(1) ) then
                     i_ne22 = 1
                     rate_ne22ag_nd = 10**rt_ne22ag(i_ne22)
                     rate_ne22an_nd = 10**rt_ne22an(i_ne22)
                  else if ( t9 >= t9_in(ndat_ne22-1) ) then
                     i_ne22 = ndat_ne22-1
                     rate_ne22ag_nd = 10**rt_ne22ag(i_ne22)
                     rate_ne22an_nd = 10**rt_ne22an(i_ne22)
                  else
                     do idt = 1, ndat_ne22 - 2
                        if ( t9_in(idt) <= t9 .and. t9 < t9_in(idt+1) ) then
                           i_ne22 = idt
                           h_intpl = ( t9 - t9_in(i_ne22) ) / ( t9_in(i_ne22 + 1) - t9_in(i_ne22) )
                           ! *** Ne22(a,g)
                           rate_ne22ag_nd = rt_ne22ag(i_ne22) &
                                   + h_intpl *( rt_ne22ag(i_ne22 + 1) - rt_ne22ag(i_ne22) ) 
                           rate_ne22ag_nd = 10.d0**rate_ne22ag_nd
                           ! *** Ne22(a,n)
                           rate_ne22an_nd = rt_ne22an(i_ne22) &
                                   + h_intpl *( rt_ne22an(i_ne22 + 1) - rt_ne22an(i_ne22) ) 
                           rate_ne22an_nd = 10.d0**rate_ne22an_nd
                           exit
                        end if
                     end do
                  end if
               end if

               !                                                             !
               !---- modification for Ne22(a,n)/(a,g) -----------------------!


               IF(T9 .GE. 0.08_r8)THEN
                  !-  NE 22 (A,G)  MG26
                  T9A   = T9/(ONE + 0.0548_r8 * T9)
                  T9A56 = T9A ** FIVE_SIXTHS
                  T9A13 = T9A ** THIRD
                  GT9   = ONE + 5._r8 * EXV(-14.791_r8/T9)
                  FT9A  = EXV(-(0.197_r8/T9A)**4.82_r8)
                  FPT9A = EXV(-(T9A/0.249_r8)**2.31_r8)
                  V(51) = (4.16D19*FPT9A/GT9*T9A56/TP32*EXV(-47.004_r8&
                    /T9A13)+2.08D16*FT9A/GT9*T9A56/TP32*EXV(-47.004_r8/T9A13))*RHO
               END IF

               ! *** new ne22ag by Joachim Gorres, M. Wiescher and G. Imbriani, 
               ! *** (+ James deBoer and Mary Beard) priv. communication 2014.
               if (ne22_michael) v(51) = rate_ne22ag_mw * rho
               ! *** Longland et al. 2012
               if (ne22_longland) v(51) = rate_ne22ag_lo * rho
               ! *** Talwar et al. 2015
               if (ne22_nd2015) v(51) = rate_ne22ag_nd * rho


               !- 1 NE 22 ( 1 ALPHA, 1 NEUT )  1 MG 25       CF 1988
               !      IF(T9 .GE. 0.08_r8)THEN
               !      T9A=T9/(ONE+0.0548_r8*T9)
               !      GTA=ONE+5.*EXV(-14.791_r8/T9)
               !      FTA=EXV(-(0.197_r8/T9A)**4.82_r8)
               !      T9A56=T9A**(5./6.)
               !      T9A13=T9A**(ONE/3.)
               !      V(52)=(4.16E+19*FTA/GTA*T9A56*TM32*EXV(-47.004_r8/T9A13)
               !     &      +1.44E-04/GTA*EXV(-5.577_r8/T9))*RHO
               !       END IF
               if ( t9 >= 0.1_r8 .and. t9 <= 10._r8 ) then
                  !     REAZIONE  NE22AN   JAEGER01 REC 0.1_r8<T9<10
                  V(52) = 4.04d0*EXP(-7.74d0/T9)&
                    +2.302d-4*T9**(-6.0d-1)*EXP(-6.14d0/T9)&
                    +6.9d+3*T9**(3.19d0)*EXP(-1.13d+1/T9)&
                    +1.881d+7*T9**(3.58d-1)*EXP(-2.67d+1/T9)
                  !     tt correction (NACRE)
                  V(52) = V(52)*(1.d0+2.674d0*EXP(-1.5025d+1/T9-3.21d-1*T9))
                  V(52)=V(52) * RHO
                  !C     REAZIONE  NE22AN   JAEGER01 UP 0.1_r8<T9<10
                  !       V(52)=3.68d0*DEXP(-7.7d0/T9)
                  !     a  +9.02d-4*T9**(-1.7d0)*DEXP(-6.31d0/T9)
                  !     a  +1.09d+4*T9**(2.853d0)*DEXP(-1.16d+1/T9)
                  !     a  +5.21d+6*T9**(1.05d0)*DEXP(-2.32d+1/T9)
                  !!     tt correction (NACRE)
                  !        V(52) = V(52)*(1.d0+2.674d0*DEXP(-1.5025d+1/T9-3.21d-1*T9))
                  !      V(52)=V(52) * RHO
                  !C     REAZIONE  NE22AN   JAEGER01 LOW 0.1_r8<T9<10
                  !       V(52)=4.55d0*DEXP(-7.781d0/T9)
                  !     a  +1.701d-10*T9**(-5.98d0)*DEXP(-6.22d0/T9)
                  !     a  +8.d+3*T9**(2.75d0)*DEXP(-1.155d+1/T9)
                  !     a  +1.003d+6*T9**(1.5d0)*DEXP(-2.3d+1/T9)
                  !!     tt correction (NACRE)
                  !        V(52) = V(52)*(1.d0+2.674d0*DEXP(-1.5025d+1/T9-3.21d-1*T9))
                  !      V(52)=V(52) * RHO
               end if

               ! *** new ne22an by Joachim Gorres, M. Wiescher and G. Imbriani, 
               ! *** (+ James deBoer and Mary Beard) priv. communication 2014.
               if (ne22_michael) v(52) = rate_ne22an_mw *rho
               if (ne22_longland) v(52) = rate_ne22an_lo *rho
               if (ne22_nd2015) v(52) = rate_ne22an_nd *rho

               !	print *,'ne22an vital',V(52)/rho
               !
               IF(T9 .GE. 0.08_r8)THEN
                  !-MG24(A,G)SI28
                  V(53)=(4.78D01/TP32*EXV(-13.506_r8/T9)+2.38D03/TP32*EXV(&
                      -15.218_r8/T9)+2.47D02*TP32*EXV(-15.147_r8/T9)+0.1_r8*(1.72D-09/&
                      TP32*EXV(-5.028_r8/T9)+1.25D-03/TP32*EXV(-7.929_r8/T9)&
                      +2.43D01/T9*EXV(-11.523_r8/T9)))/(ONE+5._r8*EXV(-15.882_r8/T9))&
                      *RHO
               END IF

               ! *** who is this?
               !*****        rates added by margaret

               !
               !------------ C 14  -------------------
               !- 1 C  14 (1 HE  4, 0 GAMMA)  1 O  18               CF1988
               V(54) = (3.375D+08/T9**2*EXV(-32.513_r8*TM13)+1.528D+09*TM23&
                                *EXV(-32.513_r8*TM13-(T9/2.662_r8)**2)*(1.0_r8+0.0128_r8*TP13-0.869_r8&
                                *TP23-0.0779_r8*T9+0.321_r8*TP43+0.0732_r8*TP53)+9.29D-08*TM32&
                         *EXV(-2.048_r8/T9)+2.77D+03*TM45*EXV(-9.876_r8/T9))* RHO

               !- 1 C  14 (1 PROT , 0 GAMMA)  1 N  15                CF1988

               V(55) = (6.80D+06/TP23*EXV(-13.741_r8/TP13-(T9/5.721_r8)**2)*(1.0_r8+0.030_r8&
                                *TP13+0.503_r8*TP23+0.107_r8*T9+0.213_r8*TP43+0.115_r8*TP53)+5.36D+03&
                             /TP32*EXV(-3.811_r8/T9)+9.82D+04/TP13*EXV(-4.739_r8/T9))&
                               * RHO * PROTOFAK

               !
               !-------------- N 15 ----------------
               !- 1 N  15 (1 ALPHA, 0 GAMMA)  1 F  19                NACRE
               V(56) = (1.10D+11/TP23*EXV(-36.214_r8/TP13-(T9/0.6_r8)**2)&
                         + 1.65D-04/TP32*EXV(-4.224_r8/T9)&
                         + 2.66_r8/TP32*EXV(-6.220_r8/T9)&
                               + 1.56D+02/TP32*EXV(-7.764_r8/T9)&
                               + 3.92D+04/TP13*EXV(-14.522_r8/T9)) * RHO


               !
               !-------------- O 17 ----------------
               !- 1 O 17 (1 ALPHA, 0 GAMMA)  1 NE 21                CF1988

               T9A = T9/(1.0_r8+0.1646_r8*T9)
               T9A13 = T9A**THIRD
               T9A56 = T9A**FIVE_SIXTHS
               GT9 = 1.0_r8+EXV(-10.106_r8/T9)/3.0_r8
               FT9A = EXV(-(0.786_r8/T9A)**3.51_r8)
               FPT9A = EXV(-(T9A/1.084_r8)**1.69_r8)

               V(57) = (1.73D+17*FPT9A/GT9*T9A56/TP32*EXV(-39.914_r8/T9A13)&
                               +3.50e+15_r8*FT9A/GT9*T9A56/TP32*EXV(-39.914_r8/T9A13))&
                               * RHO

               ! *** warning test:
               ! *** according to Descouvemont 1993, O17(a,g) reduced by three orders of magnitude
               V(57) = V(57) * 1.d-3

               !- 1 O 17 (1 ALPHA, 1 NEUT )  1 NE 20                 NACRE
               V(58) = (4.38D+17/TP23*EXV(-39.918_r8/TP13-(T9/1.1_r8)**2)&
                               + 1.73D+03/TP32*EXV(-8.55_r8/T9)&
                               + 7.50D+05*T9**1.83_r8*EXV(-13.8_r8/T9)) * RHO


               ! ---------------------------------------
               !              CARBON BURNING - Oxygen burning
               ! ---------------------------------------
               ! *** old command: if temperature is lower than 0.3_r8 GK,
               ! *** the reactions below are not calculated, and they are zero if
               ! *** considerreaction(i) is .true. in parameter.input.
               ! *** warning. FALK?


               IF(T9 .LT. T9CN) GOTO 1254 !T9CN=0.3
               T912 = TP12
               T913 = TP13
               T914 = TP14
               T915 = TP15
               T923 = TP23
               T927 = TP27
               T932 = TP32
               T935 = TP35
               T943 = TP43
               T953 = TP53
               T954 = TP54
               T972 = TP72
               TM13 = TM13
               TM23 = TM23
               TM32 = TM32
               ALPHA=0.1_r8
               !
               T9B=T9/(ONE+0.0396_r8*T9)
               T9B13=T9B**THIRD
               T9B23=T9B**TWO_THIRDS
               T9B56=T9B**FIVE_SIXTHS
               T9H=T9/(ONE+0.055_r8*T9)
               T9H13=T9H**THIRD
               T9H23=T9H**TWO_THIRDS
               T9H56=T9H**FIVE_SIXTHS
               T9X=T9/(ONE+0.067_r8*T9)
               T9X13=T9X**THIRD
               T9X23=T9X**TWO_THIRDS
               T9X56=T9X**FIVE_SIXTHS
               GT92=ONE+5*EXV(-18.960_r8/T9)
               GT94=ONE+5.*EXV(-15.882_r8/T9)
               GT96=ONE+EXV(-9.792_r8/T9)/3.+2.*EXV(-11.773_r8/T9)/3.
               GT98=ONE+1.5_r8*EXV(-5.105_r8/T9)

               !------- C12 + C12       CF 88
               C12C12=4.27D+26*T9B56/T932*EXV(-84.165_r8/T9B13-2.12D-03*T9**3)
               ! *** test: apply multiplication factor to c12c12 cf88.
               !      C12C12 = C12C12*10.d0
               ! *** modification introduced to include the new rate of Michael Wiescher.
               !     To use Michael's rate, uncomment the line indicated below (~line 1390)
               xfind    = 0._r8
               ifindpos = 0
               mask     = 0._r8
               xfac     = 0._r8

               xfind = t9
               if(xfind <= t9c12(1))then
                  C12C12W = 0._r8
               else if(xfind > t9c12(nc12t))then
                  C12C12W = 0._r8
                  print*,'WARNING -- YOU ARE ABOVE TEMPERATURE LIMIT T9C12'
                  PRINT*,'VITAL -- CHECK WHICH C12C12 YOU USE, IF YOU ARE'
                  PRINT *,'USING THE NETWORK AT THESE TEMPERATURES! T9=',XFIND
               else
                  ifindpos = minloc(t9c12(1:nc12t),mask=t9c12(1:nc12t).ge.xfind)
                  !         print *,'ifindpos',ifindpos(1)

                  xfac = (dlog10(xfind) - dlog10(t9c12(ifindpos(1)-1)))/ &
                      (dlog10(t9c12(ifindpos(1))) - dlog10(t9c12(ifindpos(1)-1)))
                  ! *** definition of ffd = rate(1:ndensdim) interpolated in temperature
                  C12C12W(1:nc12) = c12tab(1:nc12,ifindpos(1)-1) + xfac * &
                       (c12tab(1:nc12,ifindpos(1)) - &
                        c12tab(1:nc12,ifindpos(1)-1))
               end if

               ! --- Please uncomment the following line to use Michael Wiescher's new C12C12 rate
               ! *** c12c12w(1)=new recommended
               ! *** c12c12w(2)=lower limit
               ! *** c12c12w(3)=upper limit
               !      C12C12  = C12C12W(3)
               !	print *,'c12c12',c12c12,c12c12w(1:3),
               !     1  	c12c12w(2)+c12c12w(3)
               !      c12c12a =
               !      c12c12p =

               ! --- Please uncomment the following lines to use the intermediate rate between
               !     Michael Wiescher's rate above and CF88
               !      call C12C12Inter(t9, rateout)
               !      C12C12 = rateout

               ! ***    O16 + O16 (S32)  CF88         Q = 16.542_r8
               if(T9.ge.T9OO)then
                  O16O16 = 7.10d+36/T923*exp(-135.93d0/T913 &
                     -0.629d0*T923-0.445d0*T943+0.0103_r8*T9**2)
               else if (T9.lt.T9OO)then
                  O16O16 = 0.0d0
               end if
               ! ***    C12 + O16 (SI28) CF88         Q = 16.755_r8
               if(T9.ge.T9CO)then
                  C12O16 = 1.72d+31*T9H56/T932* &
                        exp(-106.594d0/T9H13)/(exp(-0.180d0*T9H**2)+ &
                        1.06d-03*dexp(2.562d0*T9H23))
               else if (T9.lt.T9CO)then
                  C12O16 = 0.0d0
               end if
               !--------- NA 23 (P,ALPHA)
               VNA23PA=(8.56D+09/T923*EXV(-20.766_r8/T913-(T9/0.131_r8)**2)&
                 *(1. +0.02_r8*T913&
                 +8.21_r8*T923+1.15_r8*T9+44.36_r8*T943+15.84_r8*T953)+4.02_r8/T932&
                 *EXV(-1.99_r8/T9)&
                 +1.18D+04/T954*EXV(-3.148_r8/T9)+8.59D+05*T943*EXV(-4.375_r8/T9)+ALPHA&
                 *3.06D-12/T932*EXV(-0.447_r8/T9))
               !--------- NA 23 (P,GAMMA)
               VNA23PG=(2.93D+08/T923*EXV(-20.776_r8/T913-(T9/0.297_r8)**2)&
                *(1._r8+.02_r8*T913+1.61_r8*T923+.226_r8*T9+4.94_r8*T943+1.76_r8*T953)&
                +9.34D+01/T932*EXV(-2.789_r8/T9)+1.89E+04/T932*EXV(-3.434_r8/T9)&
                +5.1D+04*T915*EXV(-5.510_r8/T9))/GT98

               VMGASI=(4.78D+01/T932*EXV(-13.506_r8/T9)+2.38D+3/T932*EXV(-15.218_r8/T9&
                )+2.47D+02*T932*EXV(-15.147_r8/T9)&
                                         +ALPHA*&
                (1.72D-09/T932*EXV(-5.028_r8/T9)+1.25D-03/T932*EXV(-7.929_r8/T9)&
                +2.43D+01/T9*EXV(-11.523_r8/T9)))/GT94

               VALAMG=(1.1D+08/T923*EXV(-23.261_r8/T913-(T9/0.157_r8)**2)*(ONE+0.018_r8*T913 &
                +12.85_r8*T923+1.61_r8*T9+89.87_r8*T943+28.66_r8*T953)+1.29D+02/T932*EXV(&
                -2.517_r8/T9)+5.66D+03*T972*EXV(-3.421_r8/T9)&
                              +ALPHA*&
                (2.D-05/T932*EXV(-0.858_r8/T9)+2.81D-05/T932*EXV(-0.968_r8/T9)))/GT96

               VALPSI=(1.67D+08/T923*EXV(-23.261_r8/T913-(T9/0.155_r8)**2)&
               *(ONE+0.018_r8*T913+5.81_r8*T923+0.728_r8*T9+27.31_r8*T943+8.71_r8*T953)&
               +2.20D+00/T932*EXV(-2.269_r8/T9)+1.22D+01/T932*EXV(-2.491_r8/T9)&
               +1.50D+04*T9*EXV(-4.112_r8/T9)&
                               +ALPHA*&
                 (3.34D-07/T932*EXV(-0.858_r8/T9)+5.60D-07/T932*EXV(-0.968_r8/T9)))&
                 /GT96
               !==================   PARTITION FUNCTIONS  G(T9)   =====================
               !      T92=T9**2

               !***   GNE20=ONE+EXV((-18.70_r8+1.288_r8*T9+.03468*T92)/T9)
               !***   GNE21=4.*(ONE+EXV((-3.769_r8+.1035*T9+.06498*T92)/T9))
               !***   GNE22=ONE+EXV((-13.95_r8+1.195_r8*T9+.05358*T92)/T9)
               !***   GNA23=4.0_r8*(ONE+EXV((-4.612_r8-5.623E-04*T9+.07460*T92)/T9))
               !***   GMG24=ONE+EXV((-15.14_r8+1.253_r8*T9+.04071*T92)/T9)
               !***   GMG25=6.0_r8*(ONE+EXV((-7.775_r8-.7883*T9+.1489*T92)/T9))
               !***   GMG26=ONE+EXV((-19.30_r8+.6642*T9+.1386*T92)/T9)
               !***   GAL27=6.0_r8*(ONE+EXV((-9.412_r8-.8942*T9+.1374*T92)/T9))
               !***   GSI28=ONE+EXV((-20.02_r8+1.324_r8*T9+.03104*T92)/T9)
               !**    GSI29=2.0_r8*(ONE+EXV((-14.56_r8+.3188*T9+.1298*T92)/T9))
               !***   GSI30=ONE+EXV((-23.39_r8+.3737*T9+.1576*T92)/T9)
               !***   GP31=2.0_r8*(ONE+EXV((-13.42_r8-.2143*T9+.1660*T92)/T9))

               !      IF(INETW .EQ. 1.AND. T9 .GT. T9CN)THEN
               !C- 1 C  12 ( 1 PROT, 0 GAMMA) N13 (B+) C13                    FCZ1975
               !      V(10)=(2.04D+07/TP23*EXV(-13.690_r8/TP13-(T9/1.5_r8)**2)*(ONE+0.030_r8*TP13
               !     1 +1.19_r8*TP23+0.254_r8*T9+2.06_r8*TP43+1.12_r8*TP53) + 1.08D+05/TP32
               !     2 *EXV(-4.925_r8/T9) +2.15D+05/TP32*EXV(-18.179_r8/T9))*RHO
               !      END IF

               !------------ C 12  -------------------
               !- 2 C  12 ( 0 GAMMA, 1 NEUT )  1 MG 23        CF88 + Wiescher2006, priv. comm
               ! *** BRCCN according to Dayras et al. 1977 NUC. PHYS. A 279
               BRCCN=0d0
               if(T9.lt.0.5_r8)then
                  BRCCN=0.d0
               else if(T9.ge.0.5_r8.and.T9.le.1.5_r8)then
                  BRCCN=0.859d0*dexp(-((0.766d0/T9**3)*(1.d0 + 0.0789d0*T9 + 7.74d0*T9**3)))
               else if(T9.gt.1.5_r8.and.T9.le.5.0_r8)then
                  BRCCN=0.055d0*(1.d0 - dexp(-(0.789d0*T9 - 0.976d0)))
               else if(T9.gt.5.0_r8)then
                  BRCCN=0.02d0
               end if
               !- 2 C  12 ( 0 GAMMA, 1 PROT )  1 NA 23        CF88 + Wiescher2006, priv. comm
               V(59)=(C12C12-C12C12*BRCCN) * BRCCP * RHO
               !- 2 C  12 ( 0 GAMMA, 1 ALPHA)  1 NE 20        CF88 + Wiescher2006, priv. comm
               V(60)=(C12C12-C12C12*BRCCN) * BRCCA * RHO
               !- 2 C  12 ( 0 GAMMA, 1 NEUT)  1 MG 23        CF88 + Wiescher2006, priv. comm
               V(61)=C12C12 * BRCCN * RHO

               ! *** here I include the reverse rates of c12+c12
               ! *** from public_torch. source of reverse rate: fxtsource
               v(81:83) = 0._r8
               if(t9.gt.0.1_r8)then
                  ! *** Mg23(n,C12)C12
                  v(81) = v(61) * 3.93_r8 * exp(30.16100515d0 * t9inv)
                  ! *** Na23(p,C12)C12
                  v(82) = v(59) * 3.93_r8 * exp(-25.98325915d0 * t9inv)
                  ! *** Ne20(a,C12)C12
                  v(83) = v(60) * 2.42_r8 * exp(-53.576110995d0 * t9inv)
               end if

               !-------------- O 16 ----------------
               !- 1 O 16 (0 GAMMA, 1 ALPHA)  1 C 12
               !       REV=5.13D10*T932*EXV(-83.11_r8/T9)
               !       V(23)=V(11)*REV/RHO


               !-------------- O 18 ----------------
               !- 1 O  18 ( 1 HE  4, 1 NEUT )  1 NE 21.......................CF 1988
               !     TF=T9/(ONE+0.0483_r8*T9+0.00569_r8*TP53/(ONE+0.0483_r8*T9)**(2./3.))
               !     FT=EXV(-(0.431_r8/TF)**3.89_r8)
               !     GA=ONE+5.*EXV(-23.002_r8/T9)
               !     TF13=TF**THIRD
               !     TF56=TF**FIVE_SIXTHS
               !     V8TEMP=(7.22E+17*FT/GA*TF56/TP32*EXV(-40.056_r8/TF13)
               !    & + 150.31_r8/GA*EXV(-8.101_r8/T9))*RHO
               !     V(24)=V8TEMP

               !-------------- NE 20  -------------
               !- 1 NE 20 ( 1 ALPHA, 1 PROT )  1 NA 23                         HFCZ1984
               REV=1.25_r8 * EXV(-27.606_r8 * t9inv)
               V(62)=VNA23PA * REV * RHO
               !- 1 NE 20 ( 0 GAMMA, 1 ALPHA )  1 O 16                         HFCZ1984
               !       REV=5.65D10*T932*EXV(-54.937_r8/T9)
               !       V(27)=V(12) * REV/RHO

               !---------- NE 22   ---------------------
               !- 1 NE 22 ( 1 PROT , 0 GAMMA)  1 NA 23                         CF  1987
               !        V(63)=(1.15D+09/T923*EXV(-19.475_r8/T913) + 9.77E-12/T932*EXV
               !     &  (-0.348_r8/T9) + 8.96D+03/T932*EXV(-4.34_r8/T9) + 6.52D+04/T932
               !     &  *EXV(-5.319_r8/T9) +   7.97E+05/T912*EXV(-7.418_r8/T9) +  0.1_r8
               !     &  * 1.63D-01/T932*EXV(-1.775_r8/T9))  * RHO

               !----------- NA 23  -----------------
               !- 1 NA 23 ( 1 PROT , 0 GAMMA)  1 MG 24                          FCZ1975
               !        V(65)=VNA23PG* RHO
               !- 1 NA 23 ( 1 PROT , 1 ALPHA)  1 NE 20                         HFCZ1984
               !        V(66)=VNA23PA  * RHO
               !- 1 NA 23 ( 1 ALPHA, 1 PROT )  1 MG 26                         WFHZ1978
               V(63) = (EXV(45.74_r8-(50.20_r8/T913)*(1._r8+1.343D-02*T9 &
                  +2.502D-03*T9**2-1.526D-04*T9**3)))*RHO * TM23

               !--------------- MG 24  -------------
               !- 1 MG 24 ( 1 ALPHA, 1 PROT )  1 AL 27                          FCZ1975
               REV = 1.81_r8 * EXV(-18.575_r8 * t9inv)
               V(64) = VALAMG * REV * RHO

               !- 1 MG 24 ( 1 PROT , 1 GAMMA)  1 AL 25 (B+) MG 25               CF88
               !       GT9=ONE+5.*EXV(-15.882_r8/T9)
               !       V(69)=(5.60D+8/T923*EXV(-22.019_r8/T913)*
               !     &    (ONE+0.019_r8*T913-0.173_r8*T923-0.023_r8*T9)+
               !     &    1.48D+3/T932*EXV(-2.484_r8/T9)+4.0D+3*EXV(-4.180_r8/T9))
               !     &    *  RHO/GT9

               !- 1 MG 24 ( 0 GAMMA, 1 PROT )  1 NA 23                          FCZ1975
               !       IF(T9 .GE. 0.80D0) THEN
               !       REV=7.49D10*T932*EXV(-135.679_r8/T9)
               !       V(38)=VNA23PG * REV
               !- 1 MG 24 ( 0 GAMMA, 1 ALPHA )  1 NE 20                          FCZ195
               !       REV=6.01D10*T932*EXV(-108.059_r8/T9)
               !       V(39)=V(13)*REV/(RHO)
               !       END IF

               !-------------- MG 25  -------------------
               !- 1 MG 25 ( 1 PROT , 0 GAMMA)  1 MG 26
               !        V(70)=(6.33D-09/T932*EXV(-0.668_r8/T9)+3.55D-02/T932*
               !     &        EXV(-2.23_r8/T9)+3.88D+02/T932*EXV(-3.534_r8/T9)+
               !     &        9.31D+03*EXV(-4.18_r8/T9)) *RHO

               !- 1 MG 25 ( 1 ALPHA, 1 NEUT )  1 SI 28
               IF(INETW == 1) THEN
                  TO=T9/(1._r8 + 0.063_r8*T9)
                  GO=1._r8+10._r8*EXV(-13.18_r8/T9)/3._r8
                  TO13=TO**THIRD
                  TO56=TO**FIVE_SIXTHS
                  V(65)=(3.59D+20/GO*TO56/T932*EXV(-53.41_r8/TO13))*RHO
               END IF
               !------------ MG 26  ------------------
               ! - MG 26 ( 1 PROT , 0 GAMMA)  1 AL 27                          CF 1987
               !      G5=ONE+5.*EXV(-20.99_r8/T9)
               !      V(72)=(7.39D+08/T923*EXV(-22.042_r8/T913-(T9/0.299_r8)**2)*(ONE+0.019_r8
               !     & *T913+3.61_r8*T923+0.478_r8*T9+9.78_r8*T943+3.29_r8*T953)+2.86D+03/T932
               !     &  *EXV(-3.265_r8/T9)+7.99D+04/T932*EXV(-3.781_r8/T9)
               !     &     +ALPHA*(4.33D+10/T923
               !     & *EXV(-22.042_r8/T913-(T9/0.121_r8)**2)*(ONE+0.019_r8*T913+7.31_r8*T923+0.967_r8
               !     & *T9+40.04_r8*T943+13.47_r8*T953)+3.57E+03/T932*EXV(-2.01_r8/T9)))/G5
               !     & *RHO
               !------------ AL 27  ---------------------
               !- 1 AL 27 ( 1 PROT , 0 GAMMA)  1 SI 28                         HFCZ1984
               !        V(73)=VALPSI* RHO 
               !- 1 AL 27 ( 1 PROT , 1 ALPHA)  1 MG 24                          FCZ1975
               !        V(74)=VALAMG * RHO 
               !

               !------------ O16 + O 16  -------------------
               !- 2 O  16 ( 0 GAMMA, 1 PROT )  1 P 31        CF88 + Woosley02
               V(66)=O16O16 * BROOP * RHO
               !- 2 O  16 ( 0 GAMMA, 1 ALPHA)  1 SI 28       CF88 + Woosley02
               V(67)=O16O16 * BROOA * RHO
               !- 2 O  16 ( 0 GAMMA, 1 NEUT )  1 S 31        CF88 + Woosley02
               V(68)=O16O16 * BROON * RHO

               ! *** here I include the reverse rates of c12+o16
               ! *** from public_torch. source of reverse rate: fxtsource
               ! *** P31(p,O16)O16
               v(87) = v(66) * 5.92_r8*exp(-89.0788286d0/t9)
               ! *** Si28(a,O16)O16
               v(88) = v(67) * 3.46_r8*exp(-111.3137212d0/t9)
               ! *** S31(n,O16)O16
               v(89) = v(68) * 5.92_r8 * exp(-16.8038228d0/t9)

               ! *** note: what about o16(o16,d)P30????fxtsource

               !------------ C12 + O 16  -------------------
               !- 1 O  16 ( C  12 , 1 PROT )  1 AL 27        CF88 + OLD VITAL
               V(69)=C12O16 * BRCOP * RHO
               !- 1 O  16 ( C  12 , 1 ALPHA)  1 MG 24       CF88 + OLD VITAL
               V(70)=C12O16 * BRCOA * RHO
               !- 1 O  16 ( C  12 , 1 NEUT )  1 SI 27        CF88 + OLD VITAL
               V(71)=C12O16 * BRCON * RHO
               ! *** here I include the reverse rates of c12+o16
               ! *** from public_torch. source of reverse rate: fxtsource
               ! *** Al27(p,C12)O16
               v(84) = v(69) * 1.58d0 * exp(-59.9970745d0/t9)
               ! *** Mg24(a,C12)O16
               v(85) = v(70) * 2.83d0 * exp(-78.5648345d0/t9)
               ! *** Si27(n,C12)O16
               v(86) = v(71) * 1.58d0 * exp(4.8972467d0/t9)

               !------------ NE20 + GAMMA = O 16 + alpha  ----------
               !- 1 NE 20 ( 0 GAMMA, 1 ALPHA)  1 O  16    Clayton book pag 435
               !      V(72) = 10.d0**(12.7d0 - 28.4d0/T9)   ! s^-1
               !      V(72) = V(72)
               !- 1 NE 20 ( 0 GAMMA, 1 ALPHA)  1 O  16    NACRE reverse
               v(72) = r_o16ag * 5.653d+10*TP32 * dexp(-54.886d0/T9)
               1254   CONTINUE

               if (t9 .gt. 0.05d0)then

                  !      C14(P,N)N14         CF88
                  ! *** In the cf88 website c14(p,n)n14 is given from t9 > 0.05_r8
                  !      Q68 = -0.626_r8
                  V(75) = 7.19d+05 * (1.0d0 + 0.361d0 * TP12 + 0.502d0 * T9) &
                        * dexp(-7.263d0/T9) + 3.34d+08/TP12 * dexp(-12.246d0/T9)

                  V(75) = V(75) * rho * PROTOFAK

                  !      N14(N,P)C14    CF88 reverse
                  V(76) = 7.19d+05 * (1.0d0 + 0.361d0*TP12 + 0.502d0 * T9) &
                       * 3.33d-01 + 3.34d+08/TP12 * dexp(-4.983d0/T9) &
                       * 3.33d-01
                  V(76) = V(76) * rho

               else

                  v(75) = 0.d0
                  v(76) = 0.d0

               end if

               !      C14(A,N)O17    reaclib JINA - Rauscher et al. 1994 ApJ 429, 499
               if (t9.ge.0.01d0)then
                  a0r = -9.97213e+01_r8
                  a1r = -2.09597e+01_r8
                  a2r = -2.83860e+01_r8
                  a3r =  1.68724e+02_r8
                  a4r = -2.71387e+01_r8
                  a5r =  3.43637e+00_r8
                  a6r = -4.24734e+01_r8
                  V(77) = a0r + a1r/t9 + a2r/TP13 + &
                             a3r*TP13 + a4r*t9 + a5r*TP53 + &
                             a6r*dlog(t9)

                  V(77) = dexp(V(77))*rho
               else
                  V(77) = 0.d0
               end if

               !      O17(N,A)C14         reaclib JINA - Rauscher et al. 1994 ApJ 429, 499
               if (t9 >= 0.01_r8) then
                  a0r = -1.00418e+02_r8
                  a1r =  1.30899e-01_r8
                  a2r = -2.83860e+01_r8
                  a3r =  1.68724e+02_r8
                  a4r = -2.71387e+01_r8
                  a5r =  3.43637e+00_r8
                  a6r = -4.24734e+01_r8
                  V(78) = a0r + a1r/t9 + a2r/TP13 +&
                              a3r*TP13 + a4r*t9 + a5r*TP53 +&
                              a6r*dlog(t9)
                  V(78) =  dexp(V(78)) *rho
               else
                  V(78) = 0.d0
               end if

               ! *** else connected to if statement for t9 ge 0.001_r8
               ! *** For t9 < 0.001_r8 formula are not correct.

               ! *** NOTE = IN CASE UPDATE THE 78 BELOW IN CASE VITAL IS ENLARGED!!!!
            else if (t9 < 0.001d0 ) then
               do i = 1, 78
                  if (i.ne.5)then  ! i = 5 --> be7e+
                     v(i) = 0.d0
                  end if
               end do
            end if

         end if

         ! *** be8(g,a)a from Karlsruhe nucklidekarte 2006/2007

         v(79) = 1.034548d+16  !(s^-1)


         !	do i=1,79
         !	 print *,'v(vital)',i,v(i)
         !        end do

         ! *** beta decays terrestrial to make self-consitent vital,
         ! *** despite the choice of tbetamin.

         ! *** c11(b+)b11
         v(90)    = ln2 / 1.2228d+3  !(s^-1)
         ! *** n13(b+)C13
         v(91)    = ln2 / 5.976d+2   !(s^-1)
         ! *** c14(b-)n14
         v(92)    = ln2 / 1.8070128d+11  !(s^-1)
         ! *** f17(b+)o17
         v(93)    = ln2 / 6.48d+1  !(s^-1)
         ! *** f18(b+)o18
         v(94)    = ln2 / 6.582d+3  !(s^-1)
         ! *** na22(b+)ne22
         v(95)    = ln2 / 8.2088208d+7  !(s^-1)
         ! *** mg23(b+)na23
         v(96)    = ln2 / 1.13d+1  !(s^-1)
         ! *** al26(b+)mg26 (ground state decay)
         v(97)    = ln2 / 2.2579776d+13  !(s^-1)
         ! *** si27(b+)al27
         v(98)    = ln2 / 4.16d+0  !(s^-1)
         ! *** s31(b+)p31
         v(99)    = ln2 / 2.58d+0  !(s^-1)
         ! *** o14(b+)n14
         v(100)   = ln2 / 7.059d+1  !(s^-1)
         ! *** o15(b+)n15
         v(101)   = ln2 / 1.218d+2  !(s^-1)
         ! *** na21(b+)ne21
         v(102)   = ln2 / 2.248d+1  !(s^-1)
         ! *** al25(b+)mg25
         v(103)   = ln2 / 7.18d+0  !(s^-1)
         ! *** p29(b+)si29
         v(104)   = ln2 / 4.1d+0  !(s^-1)
         ! *** p30(b+)si30
         v(105)   = ln2 / 1.5d+2  !(s^-1)

         ! *** here I include natural alpha emission to close s-process
         ! *** path after lead. For now, in ppn there is no subroutine including
         ! *** natural emission decays.
         ! *** Bi211(alpha)Tl207(beta_istantaneus)Pb207
         v(106)   = 5.383d-3    !(s^-1)
         ! *** Po210(alpha)Pb206
         v(107)   = 5.800d-8    !(s^-1)
         ! *** Po210(n,g)Pb207
         v(108)   = 4.79d+5 * rho  ! KADoNIS 31 May 2010, RT2000


         !-N(2A,G)BE9                                                NACRE REC
         HE4HE4 = 2.43D+9 * TM23 * EXP(-13.490D0 * TM13 - (T9/0.15D0)**2.D0) *&
             (1.D0 + 74.5D0 * T9) + 6.09D+5 * TM32 * DEXP(-1.054D0/T9)
         if (T9.le.0.03d0)then
            V(109) = HE4HE4 * 6.69D-12 *&
                (1.D0 - 192.D0 * T9 + 2.48D4 * T9**2 - 1.5D6 * T9**3 +&
                 4.13D7 * T9**4 - 3.90D8 * T9**5)
         else if (T9.gt.0.03d0)then
            V(109) = HE4HE4 * 2.42D-12 * (ONE - 1.52d0 * l10t9 + 0.448d0 * l10t9**2 + 0.435d0 * l10t9**3)
         end if
         V(109) = V(109) * (RHO**2.d0)

         !-Be9(g,n)2a                                                NACRE REC

         v(110) = v(109) / (rho * rho) * 5.844d+19 * t9**3 * exp(-18.258_r8 * t9inv)

         !-D(n,g)t      jina reac

         a0r = 5.365980e+00_r8
         a1r = 0.000000e+00_r8
         a2r = 0.000000e+00_r8
         a3r = 0.000000e+00_r8
         a4r = 0.000000e+00_r8
         a5r = 0.000000e+00_r8
         a6r = 7.500000e-02_r8
         V_n111 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         !resonant1
         a0r = 6.609350e+00_r8
         a1r = 0.000000e+00_r8
         a2r = 0.000000e+00_r8
         a3r = 0.000000e+00_r8
         a4r = 0.000000e+00_r8
         a5r = 0.000000e+00_r8
         a6r = 1.000000e+00_r8
         V_r111 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)


         v(111) = dexp(V_r111)*rho + dexp(V_n111)*rho


         !v(112) = D + D -> He4 + gamma

         a0r = 3.781770e+00
         a1r = 0.000000e+00_r8
         a2r = -4.261660e+00_r8
         a3r = -1.192330e-01_r8
         a4r = 7.788290e-01_r8
         a5r = -9.252030e-02_r8
         a6r = -6.666670e-01_r8
         v(112) = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(112) = dexp(v(112))*rho

         
         !v(113) = D + D -> n + He3

         a0r = 1.975000e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -4.258600e+00_r8
         a3r = 7.334690e-01_r8
         a4r = 1.718250e-01_r8
         a5r = -3.105150e-02_r8
         a6r = -6.666670e-01_r8
         v(113) = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(113) = dexp(v(113))*rho

         !v(114) = D + D -> p + H3

         a0r = 1.984700e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -4.258600e+00_r8
         a3r = 4.342430e-01_r8
         a4r = 2.386100e-01_r8
         a5r = -3.473860e-02_r8
         a6r = -6.666670e-01_r8
         v(114) = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(114) = dexp(v(114))*rho

         !v(115) = H3 + D -> He4 + n

         a0r = 2.517940e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -4.524400e+00_r8
         a3r = 3.503370e-01_r8
         a4r = 5.874700e-01_r8
         a5r = -8.849090e+00_r8
         a6r = -6.666670e-01_r8
         V_n115 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         !resonant1
         a0r = 3.934570e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -4.5244000e+00_r8
         a3r = -1.640280e+01_r8
         a4r = 1.731030e+00_r8
         a5r = -1.229660e-01_r8
         a6r = 2.313040e+00_r8
         V_r115 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(115) = dexp(V_n115)*rho + dexp(V_r115)*rho


         !v(116) = He3 + D -> He4 + p

         a0r = 4.129690e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -7.182000e+00_r8
         a3r = -1.713490e+01_r8
         a4r = 1.369080e+00_r8
         a5r = -8.144230e-02_r8
         a6r = 3.353950e+00_r8
         V_n116 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         !resonant1
         a0r = 2.468390e+01_r8
         a1r = 0.000000e+00_r8
         a2r = -7.182000e+00_r8
         a3r = 4.732880e-01_r8
         a4r = 1.468470e+00_r8
         a5r = -2.796030e+01_r8
         a6r = -6.666670e-01_r8
         V_r116 = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(116) = dexp(V_n116)*rho + dexp(V_r116)*rho

         !v(117) = pep reaction: p + p + e- --> d + n_e

         a0r = -4.364990e+01_r8
         a1r = -2.460640e-03_r8
         a2r = -2.750700e+00_r8
         a3r = -4.248770e-01_r8
         a4r =  1.598700e-02_r8
         a5r = -6.908750e-04_r8
         a6r = -2.076250e-01_r8
         v(117) = a0r + a1r/t9 + a2r/TP13 + &
                 a3r*TP13 + a4r*t9 + a5r*TP53 + &
                 a6r*dlog(t9)

         v(117) = dexp(v(117)) * rho * rho * FABE7


   end subroutine vital_calculate_rates




   subroutine vital_index(ndim)

         integer i,ndim
         integer i_pp_d,i_dp_he3,i_2he3_2phe4,i_he4he3_be7,&
               i_be7_li7,i_li7p_2he4,i_be7p_b8,i_b8_2he4,i_b8_pbe7, &
               i_be7he4_c11,i_li7he4_b11,i_b11he4_n14n,i_b11p_3he4,&
               i_c12p_n13,i_c13p_n14,i_n13p_o14,i_n14p_o15,i_n15p_he4c12,&
               i_n15p_o16,i_o16p_f17,i_o17p_he4n14,i_o17p_f18,&
               i_o18p_he4n15,i_o18p_f19,i_f19p_he4o16,i_f19p_ne20,&
               i_ne20p_na21,i_ne21p_na22,i_ne22p_na23,i_na22p_mg23,&
               i_na23p_he4ne20,i_na23p_mg24,i_mg24p_al25,i_mg25p_al26,&
               i_mg25p_mg26,i_mg26p_al27,i_al26p_si27,i_al27p_he4mg24,&
               i_al27p_si28,i_si28p_p29,i_si29p_p30,i_si30p_p31,&
               i_3he4_c12,i_c12he4_o16,i_c13he4_no16,i_n14he4_f18,&
               i_o16he4_ne20,i_o18he4_ne22,i_ne20he4_mg24,&
               i_ne21he4_nmg24,i_ne22he4_mg26,i_ne22he4_nmg25,&
               i_mg24he4_si28,i_c14he4_o18,i_c14p_n15,i_n15he4_f19,&
               i_o17he4_ne21,i_o17he4_nne20,i_2c12_pna23,&
               i_2c12_he4ne20,i_2c12_nmg23,i_ne20he4_pna23,&
               i_na23he4_pmg26,i_mg24he4_pal27,i_mg25he4_nsi28,&
               i_2o16_pp31,i_2o16_he4si28,i_2o16_ns31,i_o16c12_pal27,&
               i_o16c12_he4mg24,i_o16c12_nsi27,i_ne20_he4o16,&
               i_n13_pc12,i_f17_po16,i_c14p_nn14,i_n14n_pc14,&
               i_c14he4_no17,i_o17n_he4c14,i_be8_2he4,i_c12_3he4,&
               i_mg23n_2c12,i_na23p_2c12,i_ne20he4_2c12,i_al27p_o16c12,&
               i_mg24he4_o16c12,i_si27n_o16c12,i_p31p_2o16,&
               i_si28he4_2o16,i_s31n_2o16,i_c11_b11,i_n13_c13,&
               i_c14_n14,i_f17_o17,i_f18_o18,i_na22_ne22,i_mg23_na23,&
               i_al26_mg26,i_si27_al27,i_s31_p31,i_o14_n14,i_o15_n15,&
               i_na21_ne21,i_al25_mg25,i_p29_si29,i_p30_si30,&
               i_bi211_he4pb207,i_po210_he4pb206,i_po210n_pb207,&
               i_2anbe9,i_be9n2a,i_d_h3,i_dd_he4,i_dd_he3,i_dd_pt,&
               i_dh3_he4,i_dhe3_he4, i_pep

         common/rates_vital/i_pp_d,i_dp_he3,i_2he3_2phe4,i_he4he3_be7,&
               i_be7_li7,i_li7p_2he4,i_be7p_b8,i_b8_2he4,i_b8_pbe7, &
               i_be7he4_c11,i_li7he4_b11,i_b11he4_n14n,i_b11p_3he4,&
               i_c12p_n13,i_c13p_n14,i_n13p_o14,i_n14p_o15,i_n15p_he4c12,&
               i_n15p_o16,i_o16p_f17,i_o17p_he4n14,i_o17p_f18,&
               i_o18p_he4n15,i_o18p_f19,i_f19p_he4o16,i_f19p_ne20,&
               i_ne20p_na21,i_ne21p_na22,i_ne22p_na23,i_na22p_mg23,&
               i_na23p_he4ne20,i_na23p_mg24,i_mg24p_al25,i_mg25p_al26,&
               i_mg25p_mg26,i_mg26p_al27,i_al26p_si27,i_al27p_he4mg24,&
               i_al27p_si28,i_si28p_p29,i_si29p_p30,i_si30p_p31,&
               i_3he4_c12,i_c12he4_o16,i_c13he4_no16,i_n14he4_f18,&
               i_o16he4_ne20,i_o18he4_ne22,i_ne20he4_mg24,&
               i_ne21he4_nmg24,i_ne22he4_mg26,i_ne22he4_nmg25,&
               i_mg24he4_si28,i_c14he4_o18,i_c14p_n15,i_n15he4_f19,&
               i_o17he4_ne21,i_o17he4_nne20,i_2c12_pna23,&
               i_2c12_he4ne20,i_2c12_nmg23,i_ne20he4_pna23,&
               i_na23he4_pmg26,i_mg24he4_pal27,i_mg25he4_nsi28,&
               i_2o16_pp31,i_2o16_he4si28,i_2o16_ns31,i_o16c12_pal27,&
               i_o16c12_he4mg24,i_o16c12_nsi27,i_ne20_he4o16,&
               i_n13_pc12,i_f17_po16,i_c14p_nn14,i_n14n_pc14,&
               i_c14he4_no17,i_o17n_he4c14,i_be8_2he4,i_c12_3he4,&
               i_mg23n_2c12,i_na23p_2c12,i_ne20he4_2c12,i_al27p_o16c12,&
               i_mg24he4_o16c12,i_si27n_o16c12,i_p31p_2o16,&
               i_si28he4_2o16,i_s31n_2o16,i_c11_b11,i_n13_c13,&
               i_c14_n14,i_f17_o17,i_f18_o18,i_na22_ne22,i_mg23_na23,&
               i_al26_mg26,i_si27_al27,i_s31_p31,i_o14_n14,i_o15_n15,&
               i_na21_ne21,i_al25_mg25,i_p29_si29,i_p30_si30,&
               i_bi211_he4pb207,i_po210_he4pb206,i_po210n_pb207,&
               i_2anbe9,i_be9n2a,i_d_h3,i_dd_he4,i_dd_he3,i_dd_pt,&
               i_dh3_he4,i_dhe3_he4, i_pep


         i = 0
         i_pp_d          =   1
         i_dp_he3        =   2
         i_2he3_2phe4   =   3
         i_he4he3_be7    =   4
         i_be7_li7       =   5
         i_li7p_2he4     =   6
         i_be7p_b8       =   7
         i_b8_2he4       =   8
         i_b8_pbe7       =   9
         i_be7he4_c11    =  10
         i_li7he4_b11    =  11
         i_b11he4_n14n   =  12
         i_b11p_3he4     =  13
         i_c12p_n13      =  14
         i_c13p_n14      =  15
         i_n13p_o14       =  16
         i_n14p_o15      =  17
         i_n15p_he4c12   =  18
         i_n15p_o16      =  19
         i_o16p_f17      =  20
         i_o17p_he4n14   =  21
         i_o17p_f18      =  22
         i_o18p_he4n15   =  23
         i_o18p_f19      =  24
         i_f19p_he4o16   =  25
         i_f19p_ne20     =  26
         i_ne20p_na21    =  27
         i_ne21p_na22    =  28
         i_ne22p_na23    =  29
         i_na22p_mg23    =  30
         i_na23p_he4ne20 =  31
         i_na23p_mg24    =  32
         i_mg24p_al25    =  33
         i_mg25p_al26    =  34
         i_mg25p_mg26    =  35
         i_mg26p_al27    =  36
         i_al26p_si27    =  37
         i_al27p_he4mg24 =  38
         i_al27p_si28    =  39
         i_si28p_p29     =  40
         i_si29p_p30     =  41
         i_si30p_p31     =  42
         i_3he4_c12      =  43
         i_c12he4_o16    =  44
         i_c13he4_no16   =  45
         i_n14he4_f18    =  46
         i_o16he4_ne20   =  47
         i_o18he4_ne22   =  48
         i_ne20he4_mg24  =  49
         i_ne21he4_nmg24 =  50
         i_ne22he4_mg26  =  51
         i_ne22he4_nmg25 =  52
         i_mg24he4_si28  =  53
         i_c14he4_o18    =  54
         i_c14p_n15      =  55
         i_n15he4_f19    =  56
         i_o17he4_ne21   =  57
         i_o17he4_nne20  =  58
         i_2c12_pna23    =  59
         i_2c12_he4ne20  =  60
         i_2c12_nmg23    =  61
         i_ne20he4_pna23 =  62
         i_na23he4_pmg26 =  63
         i_mg24he4_pal27 =  64
         i_mg25he4_nsi28 =  65
         i_2o16_pp31     =  66
         i_2o16_he4si28  =  67
         i_2o16_ns31     =  68
         i_o16c12_pal27  =  69
         i_o16c12_he4mg24=  70
         i_o16c12_nsi27  =  71
         i_ne20_he4o16  =  72
         i_n13_pc12  =  73
         i_f17_po16  =  74
         i_c14p_nn14 =  75
         i_n14n_pc14 =  76
         i_c14he4_no17  =  77
         i_o17n_he4c14  =  78
         i_be8_2he4      =  79
         i_c12_3he4      =  80
         i_mg23n_2c12   =  81
         i_na23p_2c12    =  82
         i_ne20he4_2c12  =  83
         i_al27p_o16c12  =  84
         i_mg24he4_o16c12=  85
         i_si27n_o16c12  =  86
         i_p31p_2o16 =  87
         i_si28he4_2o16  =  88
         i_s31n_2o16 =  89
         i_c11_b11   =  90
         i_n13_c13   =  91
         i_c14_n14   =  92
         i_f17_o17   =  93
         i_f18_o18   =  94
         i_na22_ne22 =  95
         i_mg23_na23 =  96
         i_al26_mg26 =  97
         i_si27_al27 =  98
         i_s31_p31   =  99
         i_o14_n14   = 100
         i_o15_n15   = 101
         i_na21_ne21 = 102
         i_al25_mg25 = 103
         i_p29_si29  = 104
         i_p30_si30  = 105
         i_bi211_he4pb207= 106
         i_po210_he4pb206= 107
         i_po210n_pb207  = 108
         i_2anbe9        = 109
         i_be9n2a        = 110
         i_d_h3          = 111
         i_dd_he4        = 112
         i_dd_he3        = 113
         i_dd_pt         = 114
         i_dh3_he4       = 115
         i_dhe3_he4      = 116
         i_pep           = 117

!         if (i_be9n2a /= ndim) stop 'vital: wrong number of rates in ppn_physics.input'
         if (i_pep /= ndim) stop 'vital: wrong number of rates in ppn_physics.input'
   end subroutine vital_index


   ! *********************************************************************
   !                        --- Custom ad-hoc rates --
   ! *********************************************************************

   !  *** NOTE:
   !     This section is implemented in order to allow for quick ad-hoc rates
   !     to be applied to the code, not permanent additions or compilations.

   subroutine N14TEST(t9, rateout)
         ! *** This routine is a test routine which is an example of using an ad-hoc
         !     tabulated rate in Vital.  In the code above you should put something
         !     like:
         !     call N14TEST(t9, rateout)
         !     v(17) = rateout * RHO * PROTOFAK
         !     The rate for N14(p,g) is then set to the tabulated rate, which is
         !     interpolated for the temperature t9.
         integer :: i
         integer, parameter :: ndim = 54
         real(r8) :: t9rate(ndim), rate(ndim), lograte(ndim)
         real(r8) :: t9, rateout

         t9rate = (/0.007_r8, 0.008_r8, 0.009_r8, 0.010_r8, 0.011_r8, 0.012_r8, &
               0.013_r8, 0.014_r8, 0.015_r8, 0.016_r8, 0.018_r8, 0.020_r8, &
               0.025_r8, 0.030_r8, 0.040_r8, 0.050_r8, 0.060_r8, 0.070_r8, &
               0.080_r8, 0.090_r8, 0.100_r8, 0.110_r8, 0.120_r8, 0.130_r8, &
               0.140_r8, 0.150_r8, 0.160_r8, 0.180_r8, 0.200_r8, 0.250_r8, &
               0.300_r8, 0.350_r8, 0.400_r8, 0.450_r8, 0.500_r8, 0.600_r8, &
               0.700_r8, 0.800_r8, 0.900_r8, 1.000_r8, 1.250_r8, 1.500_r8, &
               1.750_r8, 2.000_r8, 2.500_r8, 3.000_r8, 3.500_r8, 4.000_r8, &
               5.000_r8, 6.000_r8, 7.000_r8, 8.000_r8, 9.000_r8, 10.00_r8/)
         rate = (/1.9e-26_r8, 5.6e-25_r8, 9.8e-24_r8, 1.1e-22_r8, 9.8e-22_r8, &
               6.5e-21_r8, 3.6e-20_r8, 1.7e-19_r8, 6.6e-19_r8, 2.4e-18_r8, &
               2.2e-17_r8, 1.6e-16_r8, 7.5e-15_r8, 1.4e-13_r8, 1.1e-11_r8, &
               2.2e-10_r8, 2.3e-09_r8, 1.4e-08_r8, 6.6e-08_r8, 2.4e-07_r8, &
               7.2e-07_r8, 2.18e-06_r8, 5.56e-06_r8, 1.44e-05_r8, 3.93e-05_r8, &
               0.000109_r8, 0.000292_r8, 0.00171_r8, 0.00741_r8, 0.104_r8, &
               0.584_r8, 1.94_r8, 4.67_r8, 9.04_r8, 15.1_r8, 31.4_r8, 50.5_r8, &
               70.9_r8, 90.4_r8, 108.0_r8, 145.0_r8, 174.0_r8, 202.0_r8, &
               235.0_r8, 323.0_r8, 454.0_r8, 638.0_r8, 881.0_r8, 1520.0_r8, &
               2270.0_r8, 3060.0_r8, 3800.0_r8, 4460.0_r8, 5030.0_r8/)

         !     I am interested in performing a logarithmic interpolation, so I pass
         !     the log of the rate instead of the actual rate.
         do i = 1, ndim
            lograte(i) = log10(rate(i))
         end do

         !     This function performs the interpolation and outputs the rate (rateout)
         !     for a particular temperature (t9choice)
         call rate_interp(t9rate, lograte, ndim, t9, rateout)

         rateout = 10._r8 ** rateout

         return

   end subroutine N14TEST

   ! *************************************************************************
   subroutine C12C12Inter(t9, rateout)
         ! *** This is the intermediate C12C12 rate; a geometric mean between the
         !     Michael Wiescher rate in Vital and Caughlan & Fowler 1988
         !     To use it, put into the vital subroutine code above:
         !     call C12C12Inter(t9, rateout)
         !     C12C12 = rateout

         integer ndim, i
         parameter (ndim = 49)
         double precision t9rate(ndim), rate(ndim), lograte(ndim)
         double precision t9, rateout

         t9rate = (/0.11d0, 0.121d0, 0.133d0, 0.146d0, 0.161d0, 0.177d0, &
               0.195d0, 0.214d0,0.236d0, 0.259d0, 0.285d0, 0.314d0, &
               0.345d0, 0.38d0, 0.418d0, 0.46d0, 0.505d0, 0.556d0, &
               0.612d0, 0.673d0, 0.74d0, 0.814d0, 0.895d0, 0.985d0, &
               1.08d0, 1.19d0, 1.31d0, 1.44d0, 1.59d0, 1.75d0, 1.92d0, &
               2.11d0, 2.32d0, 2.56d0, 2.81d0, 3.09d0, 3.4d0, 3.74d0, &
               4.11d0, 4.53d0, 4.98d0, 5.48d0, 6.02d0, 6.63d0, 7.29d0, &
               8.02d0, 8.82d0, 9.7d0, 10.7d0/)

         rate = (/9.81828475983d-50, 2.21439875386d-47, 4.12727078059d-45, &
               6.15573972099d-43, 8.982592702d-41, 1.0146752308d-38, &
               1.06524636834d-36, 9.85897156895d-35, 1.82590390573d-32, &
               3.48966652823d-30, 5.14442872395d-28, 5.55955173691d-26, &
               4.16263627804d-24, 2.43773419973d-22, 1.06720356735d-20, &
               3.67860381636d-19, 9.6442623973d-18, 2.15478931447d-16, &
               3.93975038259d-15, 5.9086135257d-14, 7.49916322648d-13, &
               8.19083768649d-12, 7.69680150202d-11, 6.40391449549d-10, &
               4.47377809702d-09, 3.00287819531d-08, 1.80522182076d-07, &
               1.00894968665d-06, 6.01442559034d-06, 3.76687606011d-05, &
               0.000243568732293d0,0.00155225373414d0,0.00980784746745d0, &
               0.0578178965657d0, 0.302401454093d0, 1.4927448567d0, &
               6.57911260752d0, 27.9922189182d0, 104.173075235d0, &
               364.182464175d0, 1152.38293872d0, 3162.45960711d0, &
               7876.53897039d0, 18447.2619346d0, 37229.0442071d0, &
               70462.7681714d0, 125404.571611d0, 197171.225706d0, &
               287022.015039d0/)

         !     Logarithmic interpolation: use the log of the rate.
         do i = 1, ndim
            lograte(i) = DLOG10(rate(i))
         end do

         !     This function performs the interpolation and outputs the rate (rateout)
         !     for a particular temperature (t9choice)
         call RATE_INTERP(t9rate, lograte, ndim, t9, rateout)

         rateout = 10._r8**rateout
         return
   end subroutine C12C12Inter


   ! *********************************************************************
   subroutine RATE_INTERP(t9in, ratein, ndim, t9choice, rateout)
         ! *** This subroutine performs an interpolation and returns the reaction rate
         !     for temperatures between tabulated values.  This function is to be used
         !     for ad-hoc rates in Vital in table format.
         ! *** t9in and ratein are arrays containing the temperature (in T9) and rates
         !     at those temperature.  ndim is the number of entries in the table (the
         !     logical array size).  t9choice is the temperature you wish to have the
         !     rate for.  Rateout is the final output rate.
         integer  :: ndim
         integer  :: ifindpos(1), ipos
         real(r8) :: t9in(ndim), ratein(ndim), t9choice, rateout
         real(r8) :: xfind, mask, xfac

         xfind    = 0._r8
         ifindpos = 0
         mask     = 0._r8
         xfac     = 0._r8

         if (t9choice <= t9in(1)) then
            rateout = ratein(1)
         else if (t9choice >= t9in(ndim)) then
            rateout = ratein(ndim)
         else
            ifindpos = minloc(t9in(1:ndim), mask = t9in(1:ndim) > t9choice)
            ipos = ifindpos(1)
            xfac = (t9choice - t9in(ipos-1))/(t9in(ipos) - t9in(ipos-1))
            ! *** definition of ffd = rate(1:ndensdim) interpolated in temperature
            rateout = ratein(ipos-1) + xfac*(ratein(ipos) - ratein(ipos-1))
         end if

         return
   end subroutine RATE_INTERP

   ! *********************************************************************
   subroutine CUBIC(T9,VO18)
         !..... CUBIC and related subroutines are required to interpolate
         !..... the O18(a,g) rate give in form of a Table (see BLOCK data O18gr)

         integer NG,I
         parameter(NG=53)
         integer IPOINT
         real(r8) ::  YP1,Yp2,YPN,YPN1,T9,VO18
         real(r8) ::  A(NG),B(NG),Y2(NG)
         !
         DO I=1,NG
            A(I)=TG(I)
            B(I)=DLOG10(O18AG(I))
         end do
         IPOINT=NG
         ! *** start point
         YP1 = (B(2) - B(1)) / (A(2) - A(1))
         YP2 = (B(3) - B(2)) / (A(3) - A(2))
         ! *** end point
         YPN = (B(IPOINT) - B(IPOINT-1)) / (A(IPOINT) - A(IPOINT-1))
         YPN1= (B(IPOINT-1)-B(IPOINT-2))/(A(IPOINT-1) - A(IPOINT-2))
         ! *** build spline function
         CALL SPLINE(A,B,IPOINT,YP1,YPN,Y2)
         ! *** interploate individual values
         CALL SPLINT(A,B,Y2,IPOINT,T9,VO18)
         VO18 = 10._r8**(VO18)
         !
         return

   end subroutine cubic



   function exv(x)
         real(r8) :: exv, x

         if (x >  200._r8) then
            exv = (x - 200._r8) * exp(200._r8)
         else if (x < -200._r8) then
            exv = -exp(-200._r8)/(x + 200._r8)
         else
            exv = exp(x)
         end if
   end function exv



   subroutine SPLINT(XA,YA,Y2A,N,X,Y)
         implicit double precision (A-H,O-Z)

         integer N,KLO,KHI,K
         real(r8) XA(N),YA(N),Y2A(N),X,Y,H,A,B


         KLO=1
         KHI=N
         do while (KHI-KLO.GT.1)
            K=(KHI+KLO)/2
            IF(XA(K).GT.X)THEN
               KHI=K
            ELSE
               KLO=K
            ENDIF
         end do
         H=XA(KHI)-XA(KLO)
         IF (H.EQ.0.) STOP 'Bad XA input.'
         A=(XA(KHI)-X)/H
         B=(X-XA(KLO))/H
         Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
         RETURN
   end subroutine SPLINT

   ! *********************************************************************

   subroutine SPLINE(X,Y,N,YP1,YPN,Y2)
         implicit double precision(A-H,O-Z)

         integer N,I,K
         real(r8) X(N),Y(N),Y2(N),U(N)
         real(r8) SIG,P,QN,UN,YP1,YPN

         IF (YP1.GT..99E30) THEN
            Y2(1)=0._r8
            U(1)=0._r8
         ELSE
            Y2(1)=-0.5_r8
            U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
         ENDIF
         DO I=2,N-1
            SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
            P=SIG*Y2(I-1)+2.
            Y2(I)=(SIG-1.)/P
            U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
                  /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
         end do
         IF (YPN.GT..99E30) THEN
            QN=0.
            UN=0.
         ELSE
            QN=0.5_r8
            UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
         ENDIF
         Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
         DO K=N-1,1,-1
            Y2(K)=Y2(K)*Y2(K+1)+U(K)
         end do

         return

   end subroutine spline

   ! *********************************************************************


end module vital

! *********************************************************************
! *** charged particle reaction network using formula.
! *** The list of reaction in this subr. MUST be the same as in
!     ppn_physics.input.
!
! 2 PROT    0 OOOOO  0 OOOOO   1 H   2     1.179   1
! 1 H	2     1 PROT   0 OOOOO   1 HE  3     5.494   2
! 2 HE  3   0 OOOOO  2 PROT    1 HE  4    12.860   3
! 1 HE  4   1 HE  3  0 OOOOO   1 BE  7     1.588   4
! 1 BE  7   0 OOOOO  0 OOOOO   1 LI  7     0.049   5
! 1 LI  7   1 PROT   0 OOOOO   2 HE  4    17.346   6
! 1 BE  7   1 PROT   0 OOOOO   1 B   8     0.137   7
! 1 B	8     0 OOOOO  0 OOOOO   2 HE  4    10.874   8
! 1 B	8     0 OOOOO  1 PROT    1 BE  7    -0.137   9
! 1 BE  7   1 HE  4  0 OOOOO   1 C  11     0.000  10
! 1 LI  7   1 HE  4  0 OOOOO   1 B  11     0.000  11
! 1 B  11   1 HE  4  1 NEUT    1 N  14     0.000  12
! 1 B  11   1 PROT   0 OOOOO   3 HE  4     0.000  13
! 1 C  12   1 PROT   0 OOOOO   1 N  13     3.455  14
! 1 C  13   1 PROT   0 OOOOO   1 N  14     7.551  15
! 1 N  13   1 PROT   0 OOOOO   1 N  14     4.628  16
! 1 N  14   1 PROT   0 OOOOO   1 N  15     9.059  17
! 1 N  15   1 PROT   1 HE  4   1 C  12     4.966  18
! 1 N  15   1 PROT   0 OOOOO   1 O  16    12.128  19
! 1 O  16   1 PROT   0 OOOOO   1 F  17     2.423  20
! 1 O  17   1 PROT   1 HE  4   1 N  14     1.191  21
! 1 O  17   1 PROT   0 OOOOO   1 F  18     6.867  22
! 1 O  18   1 PROT   1 HE  4   1 N  15     3.980  23
! 1 O  18   1 PROT   0 OOOOO   1 F  19     0.000  24
! 1 F  19   1 PROT   1 HE  4   1 O  16     0.000  25
! 1 F  19   1 PROT   0 OOOOO   1 NE 20     0.000  26
! 1 NE 20   1 PROT   0 OOOOO   1 NE 21     0.000  27
! 1 NE 21   1 PROT   0 OOOOO   1 NA 22     0.000  28
! 1 NE 22   1 PROT   0 OOOOO   1 NA 23     0.000  29
! 1 NA 22   1 PROT   0 OOOOO   1 NA 23     0.000  30
! 1 NA 23   1 PROT   1 HE  4   1 NE 20     0.000  31
! 1 NA 23   1 PROT   0 OOOOO   1 MG 24     0.000  32
! 1 MG 24   1 PROT   0 OOOOO   1 MG 25     0.000  33
! 1 MG 25   1 PROT   0 OOOOO   1 AL 26     0.000  34
! 1 MG 25   1 PROT   0 OOOOO   1 MG 26     0.000  35
! 1 MG 26   1 PROT   0 OOOOO   1 AL 27     0.000  36
! 1 AL 26   1 PROT   0 OOOOO   1 AL 27     0.000  37
! 1 AL 27   1 PROT   1 HE  4   1 MG 24     0.000  38
! 1 AL 27   1 PROT   0 OOOOO   1 SI 28     0.000  39
! 1 SI 28   1 PROT   0 OOOOO   1 SI 29     0.000  40
! 1 SI 29   1 PROT   0 OOOOO   1 SI 30     0.000  41
! 1 SI 30   1 PROT   0 OOOOO   1 P  31     0.000  42
! 3 HE  4   0 OOOOO  0 OOOOO   1 C  12     7.275  43
! 1 C  12   1 HE  4  0 OOOOO   1 O  16     7.162  44
! 1 C  13   1 HE  4  1 NEUT    1 O  16     0.000  45
! 1 N  14   1 HE  4  0 OOOOO   1 F  18     0.000  46
! 1 O  16   1 HE  4  0 OOOOO   1 NE 20     4.734  47
! 1 O  18   1 HE  4  0 OOOOO   1 NE 22     0.000  48
! 1 NE 20   1 HE  4  0 OOOOO   1 MG 24     9.312  49
! 1 NE 21   1 HE  4  1 NEUT    1 MG 24     2.551  50
! 1 NE 22   1 HE  4  0 OOOOO   1 MG 26    10.612  51
! 1 NE 22   1 HE  4  1 NEUT    1 MG 25    -0.480  52
! 1 MG 24   1 HE  4  0 OOOOO   1 SI 28     9.984  53
! 1 C  14   1 HE  4  0 OOOOO   1 O  18     0.000  54
! 1 C  14   1 PROT   0 OOOOO   1 N  15     0.000  55
! 1 N  15   1 HE  4  0 OOOOO   1 F  19     0.000  56
! 1 O  17   1 HE  4  0 OOOOO   1 NE 21     0.000  57
! 1 O  17   1 HE  4  1 NEUT    1 NE 20     0.000  58
! 2 C  12   0 OOOOO  1 PROT    1 NA 23     2.240  59
! 2 C  12   0 OOOOO  1 HE  4   1 NE 20     4.617  60
! 2 C  12   0 OOOOO  1 NEUT    1 MG 23     0.000  61
! 1 NE 20   1 HE  4  1 PROT    1 NA 23    -2.377  62
! 1 NA 23   1 HE  4  1 PROT    1 MG 26     1.821  63
! 1 MG 24   1 HE  4  1 PROT    1 AL 27    -1.601  64
! 1 MG 25   1 HE  4  1 NEUT    1 SI 28     2.653  65
! 2 O  16   0 OOOOO  1 PROT    1 P  31     7.677  66
! 2 O  16   0 OOOOO  1 HE  4   1 SI 28     9.593  67
! 2 O  16   0 OOOOO  1 NEUT    1 S  31     1.453  68
! 1 O  16   1 C  12  1 PROT    1 AL 27     0.000  69
! 1 O  16   1 C  12  1 HE  4   1 MG 24     0.000  70
! 1 O  16   1 C  12  1 NEUT    1 SI 27     0.000  71
! 1 NE 20   0 OOOOO  1 HE  4   1 O  16     0.000  72
! 1 N  13   0 OOOOO  1 PROT    1 C  12     0.000  73
! 1 F  17   0 OOOOO  1 PROT    1 O  16     0.000  74
! 1 C  14   1 PROT   1 NEUT    1 N  14     0.000  75
! 1 N  14   1 NEUT   1 PROT    1 C  14     0.000  76
! 1 C  14   1 HE  4  1 NEUT    1 O  17     0.000  77
! 1 O  17   1 NEUT   1 HE  4   1 C  14     0.000  78
! 1 BE  8   0 OOOOO  0 OOOOO   2 HE  4     0.000  79
! 1 C  12   0 OOOOO  0 OOOOO   3 HE  4     0.000  80
! 1 MG 23   1 NEUT   0 OOOOO   2 C  12     0.000  81
! 1 NA 23   1 PROT   0 OOOOO   2 C  12     0.000  82
! 1 NE 20   1 HE  4  0 OOOOO   2 C  12     0.000  83
! 1 AL 27   1 PROT   1 O  16   1 C  12     0.000  84
! 1 MG 24   1 HE  4  1 O  16   1 C  12     0.000  85
! 1 SI 27   1 NEUT   1 0  16   1 C  12     0.000  86
! 1 P  31   1 PROT   0 OOOOO   2 O  16     0.000  87
! 1 SI 28   1 HE  4  0 OOOOO   2 0  16     0.000  88
! 1 S  31   1 NEUT   0 OOOOO   2 O  16     0.000  89


