module slu
   use iso_c_binding
   implicit none

   interface
      subroutine slu_solve(n, nnz, isrow, shiftindex, val, colind, rowptr, rhs, ierr) bind(c)
         import :: c_int, c_double
         integer(c_int), value :: n, nnz, isrow, shiftindex
         real(c_double) :: val(*), rhs(*)
         integer(c_int) :: colind(*), rowptr(*), ierr
      end subroutine slu_solve
   end interface
end module slu
