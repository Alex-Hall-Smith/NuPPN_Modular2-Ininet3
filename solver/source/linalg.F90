module linalg
   use utils
   use solver_knobs
   use sparsekit_module, only: csrdns
   implicit none
   private
   public :: solve_Ax_equals_b_csr
contains

   !> Solve system of matrix equations Ax = b for x (A in row-compressed format)
   !> solution returned in b
   subroutine solve_Ax_equals_b_csr(Acsr,b)
#ifdef USE_SUPERLU
         use slu
#endif

#ifdef LAPACK_LIBS
#ifdef IFORT
         ! ^_^ Intel DSS MKL sparse matrix solvers
         include "mkl_dss.f77"
#endif
#endif
         type(mtx_csr) :: Acsr
         real(r8) :: b(:), perm
         real(r8), allocatable :: Adns(:,:)
         integer :: n, ierr, opt, nnz, mtx_fh, indx(nsp)
         integer, allocatable :: pvec(:)

         n = size(b)
         nnz = size(Acsr % val)

         ! print matrix to a file matrix.txt and stop
         if (.false.) then
            mtx_fh = get_unused_file_handle()
            open( mtx_fh , file = "matrix.txt" )
            call dump_mtx( n, Acsr % val, Acsr % col, Acsr % row, mtx_fh )
            close( mtx_fh )
            stop "matrix written to file matrix.txt"
         end if

         if (mat_solv_option < 4) then
            ! construct dense matrix from coordinate compressed matrix for dense solvers
            allocate(Adns(n,n))
            call csrdns( n, n, Acsr % val, Acsr % col, Acsr % row, Adns, n, ierr )
         end if

         ! matrix inversion
         if (mat_solv_option == 1) then

            ! ^_^ option 1: using Numerical Recipies 
            call ludcmp( Adns, n, n, indx, perm )
            call lubksb( Adns, n, n, indx, b )

         else if (mat_solv_option == 2) then
            !     option 2: LEQS
            !               a solver that Frank send, we tested it, and it turned
            !               out to be no improvement over option 1 for the typical
            !               application scenarios we have
            call leqs( Adns, b, n, n )

         else if (mat_solv_option == 3) then
            ! ^_^ option 3: LAPACK standard solver (from, e.g., MKL, ACML, etc.)
            call lapack_standard_solver( n, 1, Adns, n, indx, b, n, ierr )

         else if (mat_solv_option == 4) then            
            ! ^_^ option 4: Intel DSS MKL direct sparse solver.
            !     The following is a sequence of calls to the Intel DSS routines.
            !     To find a solution, one must call a series of 6 routines:
            !     dss_create
            !     dss_define_structure
            !     dss_reorder
            !     dss_factor_real
            !     dss_solve_real
            !     dss_delete
            !     The variable 'handle' is updated at each step and holds information
            !     about the particular problem you are trying to solve, so that each
            !     routine can do the appropriate procedures.

#ifndef LAPACK_LIBS
            stop 'You have selected mat_solv_option .eq. 4: Intel MKL ' // &
                  'direct sparse solver. However, LAPACK_LIBS was not ' // &
                  'defined during compilation'
#elif !defined(IFORT)
            stop 'You have selected mat_solv_option .eq. 4: Intel MKL ' // &
                  'direct sparse solver. However, you are not compiling with ' // &
                  'ARCH=INTEL.'
#else
            opt = MKL_DSS_DEFAULTS
            ierr = dss_create(handle, opt)

            opt = MKL_DSS_NON_SYMMETRIC
            ierr = dss_define_structure(handle, opt, Acsr % row, n, n, Acsr % col, nnz)

            allocate(pvec(n)) ! permutation vector
            opt = MKL_DSS_DEFAULTS
            ierr = dss_reorder(handle, opt, pvec)

            opt = MKL_DSS_INDEFINITE
            ierr = dss_factor_real(handle, opt, Acsr % val)

            opt = MKL_DSS_DEFAULTS
            allocate(rSolValues(n))
            ierr = dss_solve_real(handle, opt, b, 1, rSolValues)
            b(:) = rSolValues(:)

            opt = MKL_DSS_DEFAULTS
            ierr = dss_delete(handle,opt)
#endif

         else if ( mat_solv_option == 5 ) then
            ! ^_^ option 5: the superLU (open source) sparse solver
#ifndef USE_SUPERLU
            stop 'mat_solv_option .eq. 5 incompatible with USE_SUPERLU=NO'
#else
            call slu_solve(n, nnz, 1, 1, Acsr % val, Acsr % col, Acsr % row, b, ierr)
            !> convert back to fortran indexing
            Acsr % col(:) = Acsr % col(:) + 1
            Acsr % row(:) = Acsr % row(:) + 1

            if ( ierr /= 0 ) then
               print *, 'problem in slu_solve; ierr = ', ierr
               stop
            end if
#endif

         endif

   end subroutine solve_Ax_equals_b_csr


   !> @brief Perform a lower-upper decomposition of a matrix
   subroutine ludcmp( a, n, np, indx, d )
         use array_sizes
         integer  :: n, np, indx(n)
         real(r8) :: d, a(np,np)
         integer  :: i, imax, j, k
         real(r8) :: aamax, dum, sum, vv(n)
         real(r8), parameter :: TINY = 1.e-20_r8

         imax = 0
         d = 1._r8
         do i = 1, n
            aamax = 0._r8
            do j = 1, n
               if (dabs(a(i,j)) > aamax) aamax = dabs(a(i,j))
            end do
            if (aamax == 0._r8) stop 'singular matrix in ludcmp'
            vv(i) = 1._r8 / aamax
         end do
         do j = 1, n
            do i = 1, j-1
               sum = a(i,j)
               do k = 1, i-1
                  sum = sum - a(i,k) * a(k,j)
               end do
               a(i,j) = sum
            end do
            aamax = 0._r8
            do i = j, n
               sum = a(i,j)
               do k = 1, j-1
                  sum = sum - a(i,k) * a(k,j)
               end do
               a(i,j) = sum
               dum = vv(i) * dabs(sum)
               if (dum >= aamax) then
                  imax = i
                  aamax = dum
               endif
            end do
            if (j /= imax)then
               do k = 1, n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
               end do
               d = -d
               vv(imax) = vv(j)
            endif
            indx(j) = imax
            if (a(j,j) == 0._r8) a(j,j) = TINY
            if (j /= n) then
               dum = 1._r8 / a(j,j)
               do i = j + 1, n
                  a(i,j) = a(i,j) * dum
               end do
            endif
         end do
   end subroutine ludcmp



   !> @brief Perform a back-substitution
   subroutine lubksb(a,n,np,indx,b)
         implicit real(r8) (a-h,o-z)
         integer n,np,ii,i,j
         real(r8) :: a(np,np),b(n),sum  
         integer indx(n),ll
         ii=0
         do i=1,n
            ll=indx(i)
            sum=b(ll)
            b(ll)=b(i)
            if (ii.ne.0)then
               do j=ii,i-1
                  sum=sum-a(i,j)*b(j)
               end do
            else if (sum.ne.0.) then
               ii=i
            endif
            b(i)=sum
         end do
         do i=n,1,-1
            sum=b(i)
            if(i.lt.n)then
               do j=i+1,n
                  sum=sum-a(i,j)*b(j)
               end do
            endif
            b(i)=sum/a(i,i)
         end do
   end subroutine lubksb


   !> @brief Solves a linear system of equations a x = b via gauss jordan elimination
   !> @details a is destroyed on exit. on input b contains the right hand side, on exit it is the solution
   subroutine leqs(a,b,n,np) 
         !..declare 
         integer  ::  n, np, n1, i, j, l, imax, jj 
         real(r8) ::  a(np,np), b(np), r, c 

         !..for each column 
         n1 = n-1 
         do i=1,n 

            !..find maximum element in each row 
            r = abs(a(i,1)) 
            do j=2,n 
               c = abs(a(i,j)) 
               if (r.lt.c) r=c 
            enddo

            !..divide that row and the right hand side by the maximum element 
            !       do j=1,n 
            !        a(i,j) = a(i,j)/r 
            !       enddo

            a(i,1:n) = a(i,1:n)/r 

            b(i) = b(i)/r 
         enddo

         !..for each column, do the elimination; bail if the element is zero 
         do j=1,n1 
            l = j + 1 
            do i=l,n 
               r = -a(i,j) 
               if (r.eq.0.0_r8) goto 50 
               r = r/a(j,j) 
               !        do k=l,n 
               !         a(i,k) = a(i,k) + r*a(j,k) 
               !        enddo

               a(i,1:n) = a(i,1:n) + r*a(j,1:n)

               b(i) = b(i) + r*b(j) 
               50 continue      
            enddo
         enddo

         !..the matrix is now in triangular form; back sub 
         b(n) = b(n)/a(n,n) 
         do l=1,n1 
            i    = n-l 
            r    = 0.0d0 
            imax = i + 1 
            do j=imax,n 
               jj = i + n + 1 - j 
               r  = r + a(i,jj)*b(jj) 
            enddo
            b(i) = (b(i) - r) / a(i,i) 
         enddo
   end subroutine leqs


end module linalg
