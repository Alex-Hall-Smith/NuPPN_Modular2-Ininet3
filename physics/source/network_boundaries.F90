module network_boundaries

   !> this module deals with how the edges of the reaction network are treated. At the moment there
   !> is just one routine (natashamcloane, which is an incorrect spelling of Natascha McElhone, who
   !> is a British actress [The Truman Show ,Solaris] whose real name is actually Natasha Abigail
   !> Taylor) but there is scope to add other treatments.

   implicit none
   private

   public natashamcloane

contains

   ! find out the reaction targets, i.e. the isotope produced by reactions.  up to 4 instantaneous
   ! decays are considered after captures/photodisintegrations in the proton rich and neutron rich
   ! regions.
   ! if you want to follow the r-process or the p-process, enlarge the network!!  do not leave holes
   ! in the network (e.g.  including kr78,kr80...  without including kr79).  historically, this
   ! network is builded for the s-process.  so, check the network path to be sure to avoid these
   ! mistakes!
   subroutine natashamcloane(isotot,kk,atot,anum,ztot,znum,arrow,amin,aareac)
         use array_sizes, only: nsp
         use utils, only: r8
         implicit none
         integer :: kk, k, isotot, anum, znum, arrow, atot(nsp), ztot(nsp), areac, aareac
         real(r8) :: aanum, zznum, amin(nsp), aamin(nsp)

         aanum = dble(anum)
         zznum = dble(znum)

         areac = 0

         do k = 1, isotot
            if (atot(k) == anum .and. ztot(k) == znum) then
               arrow     = k
               areac = 1
            end if
            if (ztot(k) == znum .and. atot(k) /= 0) then
               aamin(kk) = amin(k)
            end if
         end do

         if (areac == 0) then
            do k=1,isotot
               if (atot(k) == anum .and. ztot(k) == znum-1) then
                  if (zznum/aanum > aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
               if (atot(k) == anum .and. ztot(k) == znum+1) then
                  if (zznum/aanum <= aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
            end do
         end if

         if (areac == 0) then
            do k=1,isotot
               if (atot(k) == anum .and. ztot(k) == znum-2) then
                  if ((zznum)/aanum > aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
               if (atot(k) == anum .and. ztot(k) == znum+2) then
                  if ((zznum)/aanum <= aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
            end do
         end if

         if (areac == 0) then
            do k=1,isotot
               if (atot(k) == anum .and. ztot(k) == znum-3) then
                  if ((zznum)/aanum > aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
               if (atot(k) == anum .and. ztot(k) == znum+3) then
                  if ((zznum)/aanum <= aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
            end do
         end if


         if (areac == 0) then
            do k=1,isotot
               if (atot(k) == anum .and. ztot(k) == znum-4) then
                  if ((zznum)/aanum > aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
               if (atot(k) == anum .and. ztot(k) == znum+4) then
                  if ((zznum)/aanum <= aamin(kk)) then
                     arrow     = k
                     areac = 1
                  end if
               end if
            end do
         end if

         aareac    = 0
         aamin(kk) = 0
   end subroutine natashamcloane

end module network_boundaries
