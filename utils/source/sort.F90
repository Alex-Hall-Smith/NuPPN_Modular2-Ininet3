! --------------------------------------------------------------------
! MODULE  Sorting:
!    This module can sort a set of numbers.  The method used is
! usually referred to as "selection" method.
! --------------------------------------------------------------------

MODULE  sorting
   IMPLICIT  NONE
   PRIVATE   :: FindMinimum, Swap

CONTAINS

   ! --------------------------------------------------------------------
   ! INTEGER FUNCTION  FindMinimum():
   !    This function returns the location of the minimum in the section
   ! between Start and End.
   ! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
         IMPLICIT  NONE
         double precision, DIMENSION(1:), INTENT(IN) :: x
         INTEGER, INTENT(IN)                :: Start, End
         double precision                            :: Minimum
         INTEGER                            :: Location
         INTEGER                            :: i

         Minimum  = x(Start)          ! assume the first is the min
         Location = Start             ! record its position
         DO i = Start+1, End          ! start with next elements
            IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
               Minimum  = x(i)        !      Yes, a new minimum found
               Location = i                !      record its position
            END IF
         END DO
         FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

   ! --------------------------------------------------------------------
   ! SUBROUTINE  Swap():
   !    This subroutine swaps the values of its two formal arguments.
   ! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
         IMPLICIT  NONE
         double precision, INTENT(INOUT) :: a, b
         double precision                :: Temp

         Temp = a
         a    = b
         b    = Temp
   END SUBROUTINE  Swap

   ! --------------------------------------------------------------------
   ! SUBROUTINE  Sort():
   !    This subroutine receives an array x() and sorts it into ascending
   ! order.
   ! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
         IMPLICIT  NONE
         double precision, DIMENSION(1:), INTENT(INOUT) :: x
         INTEGER, INTENT(IN)                   :: Size
         INTEGER                               :: i
         INTEGER                               :: Location

         DO i = 1, Size-1             ! except for the last
            Location = FindMinimum(x, i, Size)  ! find min from this to last
            CALL  Swap(x(i), x(Location))  ! swap this and the minimum
         END DO
   END SUBROUTINE  Sort


   subroutine indexx (n,arrin,indx)

         !     Use the Heapsort algorithm to index an array arrin of length n.
         !     Output the array indx such that arrin(indx(j)) is in ascending
         !     order for j = 1,2,...,n.  The input quantities n and arrin are
         !     not changed.

         !     This is a Numerical Recipes routine, but modified by one
         !     line to work if n equals 1.

         integer i, n, indx(n), indxt, ir, j, l
         double precision arrin(n), q

         do 11 j=1,n
            indx(j)=j
            11    continue
            if (n .eq. 1) return
            l=n/2+1
            ir=n
            10    continue
            if(l.gt.1)then
               l=l-1
               indxt=indx(l)
               q=arrin(indxt)
            else
               indxt=indx(ir)
               q=arrin(indxt)
               indx(ir)=indx(1)
               ir=ir-1
               if(ir.eq.1)then
                  indx(1)=indxt
                  return
               endif
            endif
            i=l
            j=l+l
            20      if(j.le.ir)then
               if(j.lt.ir)then
                  if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
               endif
               if(q.lt.arrin(indx(j)))then
                  indx(i)=indx(j)
                  i=j
                  j=j+j
               else
                  j=ir+1
               endif
               go to 20
            endif
            indx(i)=indxt
            go to 10
      end subroutine indexx


   END MODULE  sorting
