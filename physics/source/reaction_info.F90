module reaction_info
   use utils, only: r8, i4
   use array_sizes, only: nre
   implicit none
   private

   ! reaction types
   integer(i4), parameter, public :: &
        ! RJS 23/11/19 -- added extra reaction for beta-delayed proton emission
         num_rtypes = 22, & !< number of reaction types
         i_ng = 1, &        !< '(n,g)'
         i_gn = 2, &        !< '(g,n)'
         i_np = 3, &        !< '(n,p)'
         i_na = 4, &        !< '(n,a)'
         i_pg = 5, &        !< '(p,g)'
         i_gp = 6, &        !< '(g,p)'
         i_pn = 7, &        !< '(p,n)'
         i_pa = 8, &        !< '(p,a)'
         i_ag = 9, &        !< '(a,g)'
         i_ga = 10, &       !< '(g,a)'
         i_an = 11, &       !< '(a,n)'
         i_ap = 12, &       !< '(a,p)'
         i_bm = 13, &       !< '(-,g)' (beta minus/positron capture)
         i_ec = 14, &       !< '(+,g)' (electron capture/beta plus)
         i_bn = 15, &       !< '(b,n)' (beta-delayed neuton emission)
         i_bp = 16, &       !< '(b,p)' (beta-delated proton emission)
         ! RJS 23/11/19 -- moved all other reactions down one to add beta-delayed proton emission
         ! this is because elsewhere i > 15 (now 16) reactions are considered non-standard
         i_ba = 17, &       !< '(b,a)' (alpha-decay)
         i_ve = 18, &       !< '(v,-)' (neutrino capture)
         i_vp = 19, &       !< '(v,+)' (antineutrino capture)
         i_vsp = 20, &      !< '(v,p)' (neutrino spallation)
         i_vsn = 21, &      !< '(v,n)' (neutrino spallation)
         i_vsa = 22         !< '(v,a)' (neutrino spallation)

   ! from networksetup common block nwsetup1:
   integer(i4), public :: &
         ilabb(nre) !< reaction type label for reactions in the network (see rlabels below)

   real(r8), public :: bind_energy_diff(nre) !< reaction Q value in erg/g

   real(r8), public :: rfac(nre) !< rate multiplier factors (default to 1)

   ! reaction labels
   character(len=5), public, dimension(num_rtypes) :: rlabels !< labels for each type, e.g.  "(p,g)"


   character(len=5), public :: &
         lab(nre), & !< tag containing name for reference of a reaction rate, e.g. 'JINAC', 'KADON'
         labb(nre)   !< type of reaction, e.g, '(n,g)', '(a,p)' for each reaction in the network

         
   public reaction_info_init

contains

   subroutine reaction_info_init()
         implicit none
         integer(i4) :: i

         rlabels(i_ng) = '(n,g)'
         rlabels(i_gn) = '(g,n)'
         rlabels(i_np) = '(n,p)'
         rlabels(i_na) = '(n,a)'
         rlabels(i_pg) = '(p,g)'
         rlabels(i_gp) = '(g,p)'
         rlabels(i_pn) = '(p,n)'
         rlabels(i_pa) = '(p,a)'
         rlabels(i_ag) = '(a,g)'
         rlabels(i_ga) = '(g,a)'
         rlabels(i_an) = '(a,n)'
         rlabels(i_ap) = '(a,p)'
         rlabels(i_bm) = '(-,g)'
         rlabels(i_ec) = '(+,g)'
         rlabels(i_bn) = '(b,n)'
         ! RJS 23/11/19 -- added label for beta-delayed proton emission
         rlabels(i_bp) = '(b,p)'

         do i = 1, nre
            lab(i) = 'OOOOO'
         end do

   end subroutine reaction_info_init


end module reaction_info
