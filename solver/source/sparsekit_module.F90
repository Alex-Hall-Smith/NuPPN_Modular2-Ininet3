module sparsekit_module
   use utils, only: r8, i4

!    Selected sparsekit routines ^_^
!
!    Efstratios Gallopoulos, Youcef Saad,
!    Efficient solution of parabolic equations by Krylov approximation methods,
!    RIACS Technical Report, 90-14.
!    Noborou Kikuchi,
!    Finite element methods in mechanics,
!    Cambridge University Press, 1986.
!    David Kincaid, Thomas Oppe, John Respess, David Young,
!    ITPACKV 2C User's Guide,
!    Technical Report CNA-191.
!    Center for Numerical Analysis,
!    University of Texas at Austin, 1984.
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 3: Sorting and Searching,
!    Addison-Wesley, 1973.
!    Ole Osterby, Zahari Zlatev,
!    Direct Methods for Sparse Matrices,
!    Springer-Verlag 1983.
!    Youcef Saad,
!    Sparsekit: a basic tool kit for sparse matrix computations,
!    Technical Report, Computer Science Department,
!    University of Minnesota, June 1994.
!    Youcef Saad,
!    Analysis of some Krylov subspace approximations to the matrix exponential operator,
!    RIACS Technical Report, 90-14.
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342.
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski,
!    A testing Scheme for Subroutines Solving Large Linear Problems,
!    Computers and Chemistry,
!    Volume 5, Number 2-3, pages 91-100, 1981.
!

contains


subroutine csrcoo ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )

!*****************************************************************************80
!
!! CSRCOO converts Compressed Sparse Row to Coordinate format.
!
!  Discussion:
!
!   This routine converts a matrix that is stored in row general sparse 
!   A, JA, IA format into coordinate format AO, IR, JC. 
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) NROW, the row dimension of the matrix.
!
! job   = integer ( i4 ) serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!         (Important: note the order in the output arrays a, ja, ia. )
!         i.e., ao can be the same as a, ir can be the same as ia
!         and jc can be the same as ja.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
!
! ierr       = integer ( i4 ) error indicator.
!         ierr == 0 means normal retur
!         ierr == 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
  implicit none

  integer ( i4 ) nrow

  real ( r8 ) a(*)
  real ( r8 ) ao(*)
  integer ( i4 ) i
  integer ( i4 ) ia(nrow+1)
  integer ( i4 ) ierr
  integer ( i4 ) ir(*)
  integer ( i4 ) ja(*)
  integer ( i4 ) jc(*)
  integer ( i4 ) job
  integer ( i4 ) k
  integer ( i4 ) k1
  integer ( i4 ) k2
  integer ( i4 ) nnz
  integer ( i4 ) nzmax

  ierr = 0
  nnz = ia(nrow+1)-1

  if ( nzmax < nnz ) then
    ierr = 1
    return
  end if

  if ( 3 <= job ) then
    ao(1:nnz) = a(1:nnz)
  end if

  if ( 2 <= job ) then
    jc(1:nnz) = ja(1:nnz)
  end if
!
!  Copy backward.
!
  do i = nrow, 1, -1
    k1 = ia(i+1) - 1
    k2 = ia(i)
    do k = k1, k2, -1
      ir(k) = i
    end do
  end do

  return
end subroutine csrcoo



subroutine coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )

!*****************************************************************************80
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( i4 ) NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( i4 ) nrow

  real ( r8 ) a(*)
  real ( r8 ) ao(*)
  integer ( i4 ) i
  integer ( i4 ) iad
  integer ( i4 ) iao(nrow+1)
  integer ( i4 ) ir(*)
  integer ( i4 ) j
  integer ( i4 ) jao(*)
  integer ( i4 ) jc(*)
  integer ( i4 ) k
  integer ( i4 ) k0
  integer ( i4 ) nnz
  real ( r8 ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do
!
!  The starting position of each row.
!
  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do
!
!  Go through the structure once more.  Fill in output matrix.
!
  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do
!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1

  return
end subroutine coocsr


subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
!*****************************************************************************80
!
!! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
!
!  Discussion:
!
!    This is essentially a transposition operation.  
!
!    It is NOT an in-place algorithm.
!
!    This routine transposes a matrix stored in a, ja, ia format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) N, the order of the matrix.
!
!    Input, integer ( i4 ) JOB, indicates whether or not to fill the values of the
!    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use
!                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
!        for any other normal usage, enter ipos=1.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
!    Compressed Sparse Column format.
!
  implicit none

  integer ( i4 ) n

  real ( r8 ) a(*)
  real ( r8 ) ao(*)
  integer ( i4 ) i
  integer ( i4 ) ia(n+1)
  integer ( i4 ) iao(n+1)
  integer ( i4 ) ipos
  integer ( i4 ) j
  integer ( i4 ) ja(*)
  integer ( i4 ) jao(*)
  integer ( i4 ) job
  integer ( i4 ) k
  integer ( i4 ) next
!
!  Compute lengths of rows of A'.
!
  iao(1:n+1) = 0

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k) + 1
      iao(j) = iao(j) + 1
    end do
  end do
!
!  Compute pointers from lengths.
!
  iao(1) = ipos
  do i = 1, n
    iao(i+1) = iao(i) + iao(i+1)
  end do
!
!  Do the actual copying.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if ( job == 1 ) then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next + 1
    end do
  end do
!
!  Reshift IAO and leave.
!
  do i = n, 1, -1
    iao(i+1) = iao(i)
  end do
  iao(1) = ipos

  return
end subroutine csrcsc



subroutine csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( i4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer ( i4 ) NDNS, the dimension of the DNS array.
!
!    Output, integer ( i4 ) IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  implicit none

  integer ( i4 ) ncol
  integer ( i4 ) ndns

  real ( r8 ) a(*)
  real ( r8 ) dns(ndns,ncol)
  integer ( i4 ) i
  integer ( i4 ) ia(*)
  integer ( i4 ) ierr
  integer ( i4 ) j
  integer ( i4 ) ja(*)
  integer ( i4 ) k
  integer ( i4 ) nrow
  
  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end subroutine csrdns



subroutine dump_mtx ( n, a, ja, ia, iout )

!*****************************************************************************80
!
!! DUMP writes the matrix to a file.
!
!  Discussion:
!
!    This routine writes the matrix to a file, one row at a time in a nice 
!    readable format.  This is a simple routine which is useful for debugging.
!
!    The output unit iout will have written in it the matrix in
!    one of two possible formats (depending on the max number of
!    elements per row. the values are output with only two digits
!    of accuracy (D9.2).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer ( i4 ) IOUT, the FORTRAN output unit number.
!
  implicit none

  integer ( i4 ) n

  real ( r8 ) a(*)
  integer ( i4 ) i
  integer ( i4 ) ia(n+1)
  integer ( i4 ) iout
  integer ( i4 ) ja(*)
  integer ( i4 ) k
  integer ( i4 ) k1
  integer ( i4 ) k2
  integer ( i4 ) maxr
!
!  Select mode horizontal or vertical.
!
  maxr = 0
  do i = 1, n
    maxr = max ( maxr, ia(i+1) - ia(i) )
  end do

  if ( maxr <= 8 ) then
!
!  Able to print one row across line.
!
    do i = 1, n
      write(iout,100) i
      k1 = ia(i)
      k2 = ia(i+1) - 1
      write (iout,101) ja(k1:k2)
      write (iout,102) a(k1:k2)
    end do

  else
!
!  Unable to print one row acros line.  Do three items at a time acros line.
!
     do i = 1, n
        write(iout,200) i
        k1 = ia(i)
        k2 = ia(i+1) - 1
        write (iout,201) (ja(k),a(k), k = k1, k2)
    end do

  end if

 100  format(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  format(' col:',8(i5,6h     :))
 102  format(' val:',8(E9.2,2h :) )
 200  format(1h ,31(1h-),' row',i3,1x,31(1h-),/ &
         3('  columns :   values   *') )
 201  format(3(1h ,i5,6h    : ,D9.2,3h  *) )
  return
end subroutine dump_mtx


function getelm ( i, j, a, ja, ia, iadd, sorted )

!*****************************************************************************80
!
!! GETELM returns the element A(I,J) of a CSR matrix A. 
!
!  Discussion:
!
!    The matrix is assumed to be stored in Compressed Sparse Row (CSR) format.
!    This routine performs a binary search in the case where it is known 
!    that the elements are sorted, so that the column indices are in 
!    increasing order.  It also returns, in IADD, the address of the 
!    element A(I,J) in arrays A and JA when the search is successsful.
!    IADD is 0 if the element could not be found.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Noel Nachtigal, MIT
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) I, J, the row and column indices of the element.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(?+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ( i4 ) IADD, the address of element A(I,J) in arrays A, JA 
!    if found, zero if not found.
!
!    Input, logical SORTED, is true if the matrix is known to have its 
!    column indices sorted in increasing order.
!
!    Output, real GETELM, the value of A(I,J).
!
  implicit none

  real ( r8 ) a(*)
  real ( r8 ) getelm
  integer ( i4 ) i
  integer ( i4 ) ia(*)
  integer ( i4 ) iadd
  integer ( i4 ) ibeg
  integer ( i4 ) iend
  integer ( i4 ) imid
  integer ( i4 ) j
  integer ( i4 ) ja(*)
  integer ( i4 ) k
  logical sorted
!
!  Initialization.
!
  iadd = 0
  getelm = 0.0D+00
  ibeg = ia(i)
  iend = ia(i+1) - 1
!
!  Case where the matrix is not necessarily sorted.
!
  if ( .not. sorted ) then
!
!  Scan the row, and exit as soon as A(I,J) is found.
!
    do k = ibeg, iend
      if ( ja(k) == j ) then
        iadd = k
        go to 20
      end if
    end do

  else
!
!  Begin binary search.  Compute the middle index.
!
 10 continue

   imid = ( ibeg + iend ) / 2
!
!  Test if found.
!
     if ( ja(imid) == j ) then
        iadd = imid
        go to 20
     end if

     if ( iend <= ibeg ) then
       go to 20
     end if
!
!  else update the interval bounds.
!
     if ( j < ja(imid) ) then
        iend = imid - 1
     else
        ibeg = imid + 1
     end if

     go to 10

  end if

 20 continue

  if ( iadd /= 0 ) then
    getelm = a(iadd)
  end if

  return
end function getelm



subroutine filter ( n, job, drptol, a, ja, ia, b, jb, ib, len, ierr )

!*****************************************************************************80
!
!! FILTER copies a matrix, dropping small elements.
!
!  Discussion:
!
!    The input parameter job selects a definition of small.
!
!    This module is in place. (b,jb,ib can ne the same as
!       a, ja, ia in which case the result will be overwritten).
!
!    Contributed by David Day,  Sep 19, 1989.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) N, the row dimension of the matrix.
!
!    Input, integer ( i4 ) JOB, determines strategy chosen by caller to
!    drop elements from matrix A.
!    * 1, Elements whose absolute value is less than the drop tolerance 
!    are removed.
!    * 2, Elements whose absolute value is less than the product of the 
!    drop tolerance and the Euclidean norm of the row are removed.
!    * 3, Elements whose absolute value is less that the product of the 
!    drop tolerance and the largest element in the row are removed.
!
!    Input, real DRPTOL, the drop tolerance.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer ( i4 ) LEN, the amount of space in A and JA.
!
!    Output, real B(*), integer ( i4 ) JB(*), IB(N+1), the filtered matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ( i4 ) IERR, error flag.
!    0 indicates normal return
!    0 < IERR indicates that there is'nt enough
!    space is a and ja to store the resulting matrix.
!    IERR then contains the row number where filter stopped.
!
  implicit none

  integer ( i4 ) n

  real ( r8 ) a(*)
  real ( r8 ) b(*)
  real ( r8 ) drptol
  integer ( i4 ) ia(n+1)
  integer ( i4 ) ib(n+1)
  integer ( i4 ) ierr
  integer ( i4 ) index
  integer ( i4 ) ja(*)
  integer ( i4 ) jb(*)
  integer ( i4 ) job
  integer ( i4 ) k
  integer ( i4 ) k1
  integer ( i4 ) k2
  integer ( i4 ) len
  real ( r8 ) loctol
  real ( r8 ) norm
  integer ( i4 ) row

  index = 1

  do row = 1, n

    k1 = ia(row)
    k2 = ia(row+1) - 1
    ib(row) = index

    if ( job == 1 ) then

      norm = 1.0D+00

    else if ( job == 2 ) then

      norm = sqrt ( sum ( a(k1:k2)**2 ) )

    else

      norm = 0.0D+00

      do k = k1, k2
        if ( norm < abs ( a(k) ) ) then
          norm = abs ( a(k) )
        end if
      end do

    end if

    loctol = drptol * norm

    do k = k1, k2

      if ( loctol < abs ( a(k) ) ) then
        if ( len < index ) then
          ierr = row
          return
        end if
        b(index) = a(k)
        jb(index) = ja(k)
        index = index + 1
      end if
    end do

  end do

  ib(n+1) = index

  return
end subroutine filter


subroutine apldia ( nrow, job, a, ja, ia, diag, b, jb, ib, iw )

!*****************************************************************************80
!
!! APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in place (b, jb, ib, can be the same as
!    a, ja, ia, on entry). See comments for parameter job.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( i4 ) JOB, job indicator. Job=0 means get array b only
!    (i.e. assume that a has already been copied into array b,
!    or that algorithm is used in place. ) For all practical
!    puposes enter job=0 for an in-place call and job=1 otherwise.
!    In case there are missing diagonal elements in A,
!    then the option job =0 will be ignored, since the algorithm
!    must modify the data structure (i.e. jb, ib) in this
!    situation.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), a diagonal matrix.
!
! on return:
!
! b,
! jb,
! ib      = resulting matrix B in compressed sparse row sparse format.
!
!
! iw    = integer ( i4 ) work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer ( i4 ) nrow

  real ( r8 ) a(*)
  real ( r8 ) b(*)
  real ( r8 ) diag(nrow)
  integer ( i4 ) ia(nrow+1)
  integer ( i4 ) ib(nrow+1)
  integer ( i4 ) icount
  integer ( i4 ) ii
  integer ( i4 ) iw(*)
  integer ( i4 ) j
  integer ( i4 ) ja(*)
  integer ( i4 ) jb(*)
  integer ( i4 ) job
  integer ( i4 ) k
  integer ( i4 ) k1
  integer ( i4 ) k2
  integer ( i4 ) ko
  integer ( i4 ) nnz
  logical test
!
!  Copy integer ( i4 ) arrays into B's data structure if required.
!
  if ( job /= 0 ) then
    nnz = ia(nrow+1)-1
    jb(1:nnz) = ja(1:nnz)
    ib(1:nrow+1) = ia(1:nrow+1)
  end if
!
!  Get positions of diagonal elements in data structure.
!
  call diapos ( nrow, ja, ia, iw )
!
!  Count number of holes in diagonal and add DIAG elements to
!  valid diagonal entries.
!
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        b(iw(j)) = a(iw(j)) + diag(j)
     end if

  end do
!
!  If no diagonal elements to insert, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ib(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ib(ii)
     k2 = ib(ii+1) - 1
     ib(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = jb(k)

         if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           b(ko) = diag(ii)
           jb(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        b(ko) = a(k)
        jb(ko) = j

      end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        b(ko) = diag(ii)
        jb(ko) = ii
        iw(ii) = ko
     end if

  end do

  ib(1) = ko

end subroutine apldia 


subroutine diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
!! DIAPOS returns the positions of the diagonal elements of a sparse matrix.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) N, the row dimension of the matrix.
!
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, integer ( i4 ) IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  implicit none

  integer ( i4 ) n

  integer ( i4 ) i
  integer ( i4 ) ia(n+1)
  integer ( i4 ) idiag(n)
  integer ( i4 ) ja(*)
  integer ( i4 ) k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do i = 1, n
    do k = ia(i), ia(i+1) -1
      if ( ja(k) == i ) then
        idiag(i) = k
      end if
    end do
  end do

end subroutine diapos


subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( i4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( i4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( i4 ) n

  real ( r8 ) a(*)
  integer ( i4 ) i
  integer ( i4 ) ia(*)
  integer ( i4 ) ja(*)
  integer ( i4 ) k
  real ( r8 ) t
  real ( r8 ) x(*)
  real ( r8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do
end subroutine amux



end module sparsekit_module
