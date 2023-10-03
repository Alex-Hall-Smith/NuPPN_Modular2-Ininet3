!     A wrapper for the double-precision LAPACK solver dgesv
!
!     If LAPACK_LIBS was defined during compilation, this wrapper calls
!     dgesv.  If LAPACK_LIBS was not defined, the solver library was
!     built without LAPACK support, so dgesv is not available; the
!     wrapper explains this and stops.

subroutine lapack_standard_solver(n, nrhs, a, lda, ipiv, b, ldb, info)
      use utils, only: r8, i4
      implicit none

      integer(i4) :: n, nrhs, lda, ipiv(n), ldb, info
      real(r8) :: a(lda,lda), b(n)

#ifdef LAPACK_LIBS
      call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      if (info .ne. 0) then
         print *," WARNING nrwn: dgesv returned error info = ", info
      endif
#else
      stop "LAPACK_LIBS was not defined during compilation"
#endif
end subroutine lapack_standard_solver
