module rates
   use utils, only : r8, i4, wallclocktime, fmtprofile
   use array_sizes, only: nsp, nre
   use reaclib, only: locate_vbranch_in_reaclib_rv
   use constants, only: ONE, ZERO
#ifndef PPN
   use communication
#endif
   implicit none
   private

   ! rates array, reverse rates array, and temperature, density derivatives
   real(r8), dimension(nre), public :: v, v_rev, dvdt, dvdd

   ! reactant and product indices
   integer, dimension(nre), public  :: k1, k2, k3, k4, k5, k6, k7, k8

   ! "a, z arrays for all reactions"
   ! ant and znt are mass and charge numbers of reactant 1
   ! anf and znf are mass and charge numbers of reactant 2
   ! oh, and I have no idea about the "t" and "f"... target? f?
   integer, dimension(nre), public  :: ant, anf, znt, znf

   ! ^_^ hashing network positions
   integer, dimension(nre), public :: &
         locate_reaclib_in_vbranch, &
         locate_beta_netgen_in_vbranch, &
         locate_fuller_in_vbranch, &
         locate_oda_netgen_in_vbranch, &
         locate_lmp_netgen_in_vbranch, &
         locate_nacre_netgen_in_vbranch, &
         locate_illi_netgen_in_vbranch, &
         locate_kadonis_in_vbranch, &
         locate_jbj_in_vbranch, &
         locate_nkk_in_vbranch, &
         locate_other_in_vbranch

   public rates_init, rates_hash_locations_for_merge
#ifndef PPN
   public rates_broadcasts
#endif


contains


   subroutine rates_init()
         use nkk04, only: nkk_init
         use jbj16, only: jbj_init

         locate_reaclib_in_vbranch(:)        = -1
         locate_beta_netgen_in_vbranch(:)    = -1
         locate_fuller_in_vbranch(:)         = -1
         locate_oda_netgen_in_vbranch(:)     = -1
         locate_lmp_netgen_in_vbranch(:)     = -1
         locate_nacre_netgen_in_vbranch(:)   = -1
         locate_illi_netgen_in_vbranch(:)    = -1
         locate_kadonis_in_vbranch(:)        = -1
         locate_jbj_in_vbranch(:)            = -1
         locate_nkk_in_vbranch(:)            = -1
         locate_other_in_vbranch(:)          = -1
         locate_vbranch_in_reaclib_rv(:)     = -1

         v(:) = ZERO

         call jbj_init()
         call nkk_init()

   end subroutine rates_init



   !> @brief stores the locations of contributions from various rate sources into
   !>        vbranch array in networkII for the merge into PPN's network
   subroutine rates_hash_locations_for_merge( nrnc, ilabb )
         use reaclib, only: reaclib_ntrans
         use netgen, only: netgen_beta_ntrans, netgen_nacre_ntrans, netgen_illi_ntrans, &
               netgen_oda_ntrans, netgen_lmp_ntrans
         use kadonis, only: nkad, zkad, akad
         use fuller, only: fuller_ntrans
         use jbj16, only: jbj_ntrans
         use nkk04, only: nkk_ntrans
         use other_nuc, only: other_nuc_ntrans
         use reaction_info, only: lab, labb
         implicit none

         integer     :: kk, i, nrnc, ilabb(nre)

         do kk = 1, nrnc
            select case(lab(kk))
            case('BASEL', 'JINAR', 'JINAC', 'JINAV')
               i = reaclib_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               ! note: I just consider the case with istate = 1 in reaclib_ntrans for now.
               locate_reaclib_in_vbranch(kk) = i
               locate_vbranch_in_reaclib_rv(i) = kk
            case('NETB1')
               i = netgen_beta_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_beta_netgen_in_vbranch(kk) = i
            case('FFW85')
               i = fuller_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_fuller_in_vbranch(kk) = i
            case('ODA94')
               i = netgen_oda_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_oda_netgen_in_vbranch(kk) = i
            case('LMP00')
               i = netgen_lmp_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_lmp_netgen_in_vbranch(kk) = i
            case('NACRL','NACRR','NACRU')
               i = netgen_nacre_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_nacre_netgen_in_vbranch(kk) = i
            case('ILI01')
               i = netgen_illi_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_illi_netgen_in_vbranch(kk) = i
            case('KADON')
               do i = 1, nkad
                  if (labb(kk) == '(n,g)' .and. znt(kk) == zkad(i).and. ant(kk) == akad(i)) then
                     locate_kadonis_in_vbranch(kk) = i
                  end if
               end do
            case('JBJ16')
               i = jbj_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_jbj_in_vbranch(kk) = i
            case('NKK04')
               i = nkk_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_nkk_in_vbranch(kk) = i
            case('OOOOO','VITAL','MPG06','RVRSE','NTRNO')
               continue  ! blank, vital, a-decay neutrino or reverse reaction
            case default
               ! from other_nuc module
               i = other_nuc_ntrans(ant(kk),znt(kk),1,ilabb(kk))
               locate_other_in_vbranch(kk) = i
            end select
         end do

   end subroutine rates_hash_locations_for_merge



#ifndef PPN
   subroutine rates_broadcasts()
         return

         call broadcast(locate_reaclib_in_vbranch)
         call broadcast(locate_beta_netgen_in_vbranch)
         call broadcast(locate_fuller_in_vbranch)
         call broadcast(locate_oda_netgen_in_vbranch)
         call broadcast(locate_lmp_netgen_in_vbranch)
         call broadcast(locate_nacre_netgen_in_vbranch)
         call broadcast(locate_illi_netgen_in_vbranch)
         call broadcast(locate_kadonis_in_vbranch)
         call broadcast(locate_jbj_in_vbranch)
         call broadcast(locate_nkk_in_vbranch)
         call broadcast(locate_other_in_vbranch)
         call broadcast(locate_vbranch_in_reaclib_rv)
         call broadcast(k1)
         call broadcast(k2)
         call broadcast(k3)
         call broadcast(k4)
         call broadcast(k5)
         call broadcast(k6)
         call broadcast(k7)
         call broadcast(k8)
         call broadcast(ant)
         call broadcast(znt)
         call broadcast(znf)

   end subroutine rates_broadcasts
#endif




end module rates
