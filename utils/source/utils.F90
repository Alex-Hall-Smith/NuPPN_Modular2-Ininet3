!> @brief Utilities for the physics and solver packages and frames
!
!> details This module provides core routines and variables that are used in many of the physics and solver packages and the frames.
!> For example, definitions of data kinds, timing routines, hashing and lookup functions and so on.
!
!> @todo interpolation should go into its own module
module utils
   use array_sizes
   use sorting, only: Sort, indexx
   implicit none

   ! endianness check
   logical, parameter :: bigend = ichar(transfer(1,'a')) == 0 !< true for big endian systems

   ! ^_^ kinds a la IEEE754 via Reinecke
   integer, parameter :: &
         i4 = selected_int_kind  (9), &      !< four byte signed integer kind
         i8 = selected_int_kind  (15), &      !< eight byte signed integer kind
         r4 = selected_real_kind ( 6, 37), & !< four byte real kind (ieee754)
         r8 = selected_real_kind (15,307)    !< eight byte real kind (ieee754)

   ! ^_^ file handles
   integer, parameter :: &
         table_input = 15, &
         summary_output = 124, &
         restart_log = 235

   integer :: &
         table_output

   ! ^_^ integer constants (should probably move to physics)
   integer, parameter :: &
         irate = 1, & !< index of rate in rate arrays
         idrdt = 2, & !< index of derivative of rate wrt temperature in rate arrays
         idrdd = 3    !< index of derivative of rate wrt density in rate arrays

   real(r8) :: ispetime

   type :: mtx_csr !< CSR matrix
      real(r8), allocatable :: val(:) !< non-zero matrix elements
      integer, allocatable  :: col(:) !< columns
      integer, allocatable  :: row(:) !< rows
   end type mtx_csr

   type :: mtx_csc !< CSC matrix
      real(r8), allocatable :: val(:) !< non-zero matrix elements
      integer, allocatable  :: col(:) !< columns
      integer, allocatable  :: row(:) !< rows
   end type mtx_csc

   type :: mtx_coo !< COO matrix
      real(r8), allocatable :: val(:) !< non-zero matrix elements
      integer, allocatable  :: col(:) !< columns
      integer, allocatable  :: row(:) !< rows
   end type mtx_coo

   type :: int_hash_table !< integer hash table
      integer              :: nbuckets      !< length of hash table
      integer              :: nvalsinbucket !< number of collisions in hash table
      integer, allocatable :: values(:,:)   !< hash table
   end type int_hash_table

   type :: char_hash_table !< character hash table
      integer              :: nbuckets        !< length of hash table
      integer              :: nvalsinbucket   !< number of collisions in hash table
      character(len=5), allocatable :: values(:,:) !< hash table
   end type char_hash_table

   type :: netgen_data !< data type used for storing data from the NETGEN repository
      ! ^_^ rate data
      real(r8) :: &
            tgrid(ngrid), &               !< temperature grid
            l10t8grid(ngrid), &           !< log10 temperature grid
            aNegrid(nre1, 11), &          !< electron density grid
            vgrid(ngrid, nre1, 11), &     !< rates
            l10vgrid(ngrid, nre1, 11), &  !< log10(rates)
            flag(nre1), &                 !< -ve for T-only dependent reactions; |flag| indicates reaction type
            Qrad(nre1), &                 !< Q value
            Qnu(nre1)                     !< average neutrino energy (for weak reactions)?
      real(r8), allocatable :: rates(:)   !< interpolated rates to return to ppn
      integer  :: kgrid             !< number of temperatures for which there is actually data
      integer  :: ireac             !< no idea
      ! ^_^ reactants                        
      real(r8), dimension(nre1) :: &
            a1, a2, a3, a4                !< mass numbers
      integer, dimension(nre1)  :: &
            n1, n2, n3, n4,  &            !< stoichiometric factors
            ndata, &                      !< how many data points in electron density (I think...)
            zt, zf, at, af                !< initial, final z and a (for merging in ppn?)
      character(len=6), dimension(nre1) :: &
            z1, z2, z3, z4                !< isotope names
      character(len=2), dimension(nre1) :: &
            elt, elf                      !< initial, final element names (for merging in ppn?)
      character(len=5) :: ref             !< reference for these rates
   end type netgen_data

   integer, parameter      :: bigprime = 5381 !< a big prime number for hashing

   ! ^_^ timing stuff
   character(len=30) :: fmtprofile = "(A50,ES20.8)"

   !> Create a hash table
   interface create_hash_table
      module procedure  crt_1d_1d_char_ht, crt_2d_1d_int_ht
   end interface

   !> Deallocate a hash table
   interface free_hash_table
      module procedure  free_int_hash_table, free_char_hash_table
   end interface

   !> Generic interface to check sum of isotopic mass fractions does not deviate from 1 by more than some tolerance
   interface check_mass_conservation
      module procedure check_1d_mass_cons, check_2d_mass_cons
   end interface


contains

   !> allocate CSR matrix
   subroutine allocate_csr_matrix( mtx, nnz, nrow )
         type(mtx_csr) :: mtx
         integer :: nnz, nrow
         allocate( &
               mtx % col(nnz), &
               mtx % val(nnz), &
               mtx % row(nrow+1) &
               )
   end subroutine allocate_csr_matrix

   !> allocate COO matrix
   subroutine allocate_coo_matrix( mtx, nnz )
         type(mtx_coo) :: mtx
         integer :: nnz
         allocate( &
               mtx % col(nnz), &
               mtx % val(nnz), &
               mtx % row(nnz) &
               )
   end subroutine allocate_coo_matrix

   !> free csc matrix
   subroutine free_csc_matrix( mtx )
         type(mtx_csc) :: mtx
         deallocate(mtx % val, mtx % col, mtx % row)
   end subroutine free_csc_matrix

   !> free csr matrix
   subroutine free_csr_matrix( mtx )
         type(mtx_csr) :: mtx
         deallocate(mtx % val, mtx % col, mtx % row)
   end subroutine free_csr_matrix


   !> free coo matrix
   subroutine free_coo_matrix( mtx )
         type(mtx_coo) :: mtx
         deallocate(mtx % val, mtx % col, mtx % row)
   end subroutine free_coo_matrix

   !---
   ! DESCRIPTION:
   !> @brief Get system time in double precision
   !> @return wallclocktime (system time in double precision)
   !---
   function wallclocktime()
         real(r8) :: wallclocktime
         integer(kind=8), save :: firstcount = -huge(0)
         integer(kind=8) :: cnt, count_rate, count_max

         call system_clock(cnt, count_rate, count_max)

         if (firstcount == -huge(0)) firstcount = cnt
         if (cnt < firstcount) firstcount = firstcount - count_max
         wallclocktime = (dble(cnt - firstcount)) / count_rate
   end function wallclocktime


   !---
   ! DESCRIPTION:
   !> @brief Creates hash table for a 1d character array
   !> @detail A hash table is essentially a key--value map, like a dictionary in python
   !>     or (I believe) C++. Given a string, say, '  he4', it will provide fast
   !>     lookup of that string in a character array by using the hash function
   !>     (below) to convert the string into an (ideally unique) integer, and then
   !>     compress that integer into a range (0,...,N) by doing mod(integer,N) (in
   !>     our case, mod(integer, N) + 1 because of Fortran indexing). Naturally,
   !>     it's hard to get a perfect hash function so that all of the keys have a
   !>     *unique* hash value, and there are sometimes 2 keys with the same
   !>     compressed hash value. This is called a collision. A good hash function
   !>     minimizes collisions, but some are inevitable. So, sometimes a hashed key
   !>     (e.g. hash('  he4')) gives an index, e.g., 5, that is also given by, say,
   !>     hash('  u235'). In the case of such a collision, both values are stored
   !>     in the hash array under bucket number 5. Then, at the time we want to get
   !>     the index of '  he4', we lookup in the hash table and get 5, and then we
   !>     only have to compare '  he4' with '  he4' and '  u235' in order to get
   !>     the real index. This is at best 1 hash and 1 lookup, and for a decent
   !>     hash function will never get close to O(n).
   !> @return hash table
   !---
   function crt_1d_1d_char_ht(array_to_hash)
         type(int_hash_table) :: crt_1d_1d_char_ht
         character(len=5)     :: array_to_hash(:)
         integer              :: ncollisions(bigprime+1)
         integer              :: max_values_in_bucket
         integer              :: i, j, hashed
         integer              :: n

         n = size(array_to_hash)
         ncollisions(:) = 0
         allocate(crt_1d_1d_char_ht% values(bigprime+1,1))
         crt_1d_1d_char_ht% values(:,:) = -1
         ! ^_^ count max number of collisions
         do i = 1, n
            hashed = hash(array_to_hash(i))
            if (crt_1d_1d_char_ht% values(hashed,1) /= -1) then
               ncollisions(hashed) = ncollisions(hashed) + 1
            end if
            crt_1d_1d_char_ht% values(hashed,1) = i
         end do

         max_values_in_bucket = maxval(ncollisions) + 1
         deallocate(crt_1d_1d_char_ht% values)
         allocate(crt_1d_1d_char_ht% values(bigprime+1,max_values_in_bucket))

         ! ^_^ write hash array
         crt_1d_1d_char_ht% values(:,:) = -1
         do i = 1, n
            hashed = hash(array_to_hash(i))
            do j = 1, ncollisions(hashed) + 1
               if (crt_1d_1d_char_ht% values(hashed,j) == -1) exit
            end do
            crt_1d_1d_char_ht% values(hashed,j) = i
         end do

         crt_1d_1d_char_ht% nbuckets      = bigprime+1
         crt_1d_1d_char_ht% nvalsinbucket = max_values_in_bucket

   end function crt_1d_1d_char_ht


   !---
   ! DESCRIPTION:
   ! 
   !> @brief DJB2 Hashing function
   !> @details This function takes a variable-length string and turns it into an integer using the djb2 hash function, or something
   !> close to it anyway.
   !
   !> @returns integer (hashed string)
   !---
   integer function hash(str)
         character(len=*),intent(in) :: str !< string to be hashed
         integer :: i, j, magic_numb
         data magic_numb / Z'5D7A9F43' /

         hash = 0

         do i = 1, len(str)
            j = mod(i-1, 4) * 8
            hash = ieor( hash, ishft( ichar( str(i:i) ), j ) ) 
            hash = abs( ieor( hash, magic_numb ) )
            hash = mod(hash, bigprime) + 1
         end do
   end function hash


   !---
   ! DESCRIPTION:
   !
   !> @brief Provide a unique integer for a pair of integers
   ! 
   !> @details provide a unique integer for a pair of integers, i, j; useful for compressing sparse 2d arrays into 1d arrays.
   !  array1_to_hash and array2_to_hash are the arrays of the i and j values, respectively and must be of the same length. Assign
   !  the unique integer to the corresponding entry in the values array
   !---
   function crt_2d_1d_int_ht(array1_to_hash, array2_to_hash, values)
         type(char_hash_table) :: crt_2d_1d_int_ht
         integer          :: array1_to_hash(:), array2_to_hash(:)
         integer          :: ncollisions(bigprime+1)
         integer          :: max_values_in_bucket
         integer          :: i, j, hashed
         integer          :: n
         character(len=5) :: values(:)

         n = size(array1_to_hash)
         ncollisions(:) = 0
         allocate(crt_2d_1d_int_ht% values(bigprime+1,1))
         crt_2d_1d_int_ht% values(:,:) = 'XXXXX'
         ! ^_^ count max number of collisions
         do i = 1, n
            hashed = Cantor_hash(array1_to_hash(i), array2_to_hash(i))
            if (crt_2d_1d_int_ht% values(hashed,1) /= 'XXXXX') then
               if (array1_to_hash(i) /= 0 .and. array2_to_hash(i) /= 0) then
                  ncollisions(hashed) = ncollisions(hashed) + 1
               end if
            end if
            crt_2d_1d_int_ht% values(hashed,1) = values(i)
         end do

         max_values_in_bucket = maxval(ncollisions) + 1
         deallocate(crt_2d_1d_int_ht% values)
         allocate(crt_2d_1d_int_ht% values(bigprime+1,max_values_in_bucket))

         ! ^_^ write hash array
         crt_2d_1d_int_ht% values(:,:) = 'XXXXX'
         do i = 1, n
            hashed = Cantor_hash(array1_to_hash(i), array2_to_hash(i))
            do j = 1, ncollisions(hashed) + 1
               if (crt_2d_1d_int_ht% values(hashed,j) == 'XXXXX') exit
            end do
            crt_2d_1d_int_ht% values(hashed,j) = values(i)
         end do


         crt_2d_1d_int_ht% nbuckets      = bigprime+1
         crt_2d_1d_int_ht% nvalsinbucket = max_values_in_bucket

   end function crt_2d_1d_int_ht


   !---
   ! DESCRIPTION:
   !> @brief Cantor pairing function
   !
   !> @details create a hash integer from two integers using the Cantor pairing funtion and compression to bigprime
   !
   !> @returns integer hash of integer pair (a,b)
   !---
   integer function Cantor_hash(a, b)
         integer :: &
               a, & !< first integer of pair to be hashed
               b    !< second integer of pair to be hashed

         Cantor_hash = (a + b) * (a + b + 1) / 2 + b

         ! ^_^ compress
         Cantor_hash = mod(Cantor_hash, bigprime) + 1
   end function Cantor_hash


   !---
   ! DESCRIPTION:
   !> @brief deallocate a hash table (instance of the \p int_hash_table type
   !---
   subroutine free_int_hash_table(table)
         type(int_hash_table) :: table !< table hash table (instance of \p int_hash_table) to be freed

         if (allocated(table% values)) then

            deallocate(table% values)
            table% nbuckets      = -1
            table% nvalsinbucket = -1
         else
            stop 'error: tring to free an already empty hash table'
         end if

   end subroutine free_int_hash_table


   !---
   ! DESCRIPTION:
   !> @brief Deallocate a hash table (instance of the \p char_hash_table type
   !---
   subroutine free_char_hash_table(table)
         type(char_hash_table) :: table !< table hash table (instance of \p char_hash_table) to be freed

         if (allocated(table% values)) then

            deallocate(table% values)
            table% nbuckets      = -1
            table% nvalsinbucket = -1
         else
            stop 'error: tring to free an already empty hash table'
         end if
   end subroutine free_char_hash_table


   !---
   ! DESCRIPTION:
   !> @brief Calculate deviation from 1 of the mass fraction sum
   !> @details i.e. \f$ 1-\sum_i X_i \f$
   !
   !> @param[in] sort whether or not to sort the abundances before summing (more accurate but much slower...)
   !> @param[in] yps array containing isotopic mass fractions
   !> @param[out] residual \f$ 1-\sum_i X_i \f$
   !---
   subroutine check_1d_mass_cons(sort, yps, residual)

         real(r8) :: yps(:), residual
         logical :: sort

         if (sort) then
            residual = abs(1._r8 - sortedsum(yps))
         else
            residual = abs(1._r8 - sum(yps))
         end if
   end subroutine check_1d_mass_cons


   !---
   ! DESCRIPTION:
   !> @brief Calculate deviation from 1 of the mass fraction sum for all mass zones
   !> @details i.e. \f$ 1-\sum_i X_i \f$
   !
   !> @param[in] sort whether or not to sort the abundances before summing (more accurate but much slower...)
   !> @param[in] yps array containing isotopic mass fractions for every mass zone
   !> @param[out] residual \f$ \max(1-\sum_i X_i) \f$
   !---
   subroutine check_2d_mass_cons(sort, yps, m, residual)
         real(r8) :: yps(:,:), residual
         integer :: i, m
         logical :: sort

         residual = 0._r8
         do i = 1, m
            if (sort) then
               residual = max(residual, abs(1._r8 - sortedsum(yps(i,:))))
            else
               residual = max(residual, abs(1._r8 - sum(yps(i,:))))
            end if
         end do
   end subroutine check_2d_mass_cons


   !---
   ! DESCRIPTION:
   !> @brief Renormalise abundances to ensure \f$ \sum_i X_i = 1 \f$
   !
   !> @param[in,out] yps array containing isotopic mass fractions
   !> @param[in] residual \f$ \max(1-\sum_i X_i) \f$
   !> @param[in] message information about why the sum of isotopic mass fractions may have deviated from 1
   !---
   subroutine renorm_1d(yps, residual, message)
         real(r8) :: yps(:)
         real(r8) :: residual
         character(len=*) :: message

         if (residual > 1e-4_r8) then
            write(*,fmtprofile) 'WARNING: mass cons. residual after ' // trim(message) // ': ', residual
            write(*,*) '***RENORMALISING***'
         end if

         yps(:) = yps(:) / ksum(yps)
   end subroutine renorm_1d


   !---
   ! DESCRIPTION:
   !> @brief Renormalise abundances to ensure \f$ \sum_i X_i = 1 \f$
   !
   !> @param[in,out] yps array containing isotopic mass fractions for each mass zone
   !> @param[in] residual \f$ \max(1-\sum_i X_i) \f$
   !> @param[in] message information about why the sum of isotopic mass fractions may have deviated from 1
   !---
   subroutine renorm_2d(yps, m, max_residual, message)
         real(r8) :: yps(:,:)
         real(r8) :: max_residual
         integer :: i, m
         character(len=*) :: message

         if (message /= "") then
            write(*,fmtprofile) 'WARNING: mass cons. max. residual after ' // trim(message) // ': ', max_residual
            write(*,*) '***RENORMALISING***'
         end if

         do i = 1, m
            yps(i,:) = yps(i,:) / ksum(yps(i,:))
         end do
   end subroutine renorm_2d


   !---
   ! DESCRIPTION:
   !> @brief Sort isotopic mass fractions and calculate \f$ \sum_i X_i \f$, for better accuracy (minimizing truncation errors)
   !
   !> @details In reality, this is really inefficient and should be improved before usage
   !
   !> @param[in] array array to sort and calculate sum of
   !> @return sortedsum sum of sorted array
   !---
   function sortedsum(array, mask)
         real(r8) :: sortedsum, array(:)
         integer, allocatable :: indx(:)
         integer :: n, i
         logical, optional :: mask(:)

         n = size(array)
         allocate(indx(n))
         call indexx(n, array, indx)
         sortedsum = 0._r8
         do i = 1, n
            if ( present(mask) ) then
               if ( mask(i)) then
                  sortedsum = sortedsum + array(indx(i))
               end if
            else
               sortedsum = sortedsum + array(indx(i))
            end if
         end do

   end function sortedsum

   !> calculate abar and zbar, the average ion weight and ion charge
   subroutine calculate_azbar(x, a, z, abar, zbar, mask)
         use array_sizes, only: nsp
         implicit none
         logical :: mask(nsp)       !< .true. where isotope should be considered
         real(r8) :: &
               x(nsp), &            !< isotopic mass fractions
               a(nsp), &            !< mass numbers
               z(nsp), &            !< proton numbers
               abar, &              !< average ion weight
               zbar, &              !< average ion charge
               y(nsp), &                 !< mole fractions
               ytot                 !< total mole fraction for averaging
         integer(i4) :: i

         where(mask)
            y = x / a
         end where

         ytot = sum(y, mask = mask)
         abar = sum(a*y, mask = mask) / ytot
         zbar = sum(z*y, mask = mask) / ytot
   end subroutine calculate_azbar

   !---
   ! DESCRIPTION:
   !
   !> @brief Calculate electron fraction Ye from the neutron excess
   !
   !> @details The neutron excess is calculated by \f$ \eta = \sum_i (A_i - 2Z_i) \frac{X_i}{A_i}  \f$
   !> from which the electron fraction is calculated: \f$ Y_\mathrm{e} = \frac{1-\eta}{2} \f$
   !> This method assumes mass conservation, therefore it avoids erroneous \f$ Y_\mathrm{e} \f$ values when there are deviations
   !> from \f$ \sum_i X_i = 1 \f$
   !---
   subroutine calculate_ye(x, a, z, ye, mask)
         use array_sizes, only: nsp
         integer :: i
         real(r8) :: &
               x(nsp), &            !< isotopic mass fractions
               a(nsp), &            !< mass numbers
               z(nsp), &            !< proton numbers
               eta(nsp), &          !< neutron excess
               ye                   !< electron fraction
         logical :: mask(nsp)       !< .true. where isotope should be considered in the calculation of \f$ Y_\mathrm{e} \f$

         eta(:) = 0._r8
         where(mask)
            eta = (a - 2._r8 * z) * x / a
         end where
         ye = (1._r8 - sortedsum(eta)) * 0.5_r8 ! This is the effective Ye
   end subroutine calculate_ye



   subroutine check_file_exists(filename)
         character(len=*) :: filename
         logical :: file_exists
         inquire( file = trim(filename), exist = file_exists )
         if (.not. file_exists) then
            write(*,*) "file does not exist: ", filename
            stop "error: missing file"
         end if
   end subroutine check_file_exists


   !---
   ! DESCRIPTION:
   !
   !> @brief Convert string to all-caps
   !
   !> @details Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012); Original author: Clive Page
   !
   !> @returns capitalized string
   !---
   function to_upper(strIn) result(strOut)
         character(len=*), intent(in) :: strIn !< string to be capitatized
         character(len=len(strIn)) :: strOut
         integer :: i,j

         do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
            else
               strOut(i:i) = strIn(i:i)
            end if
         end do
   end function to_upper


   !---
   ! DESCRIPTION:
   !
   !> @brief return a file handle that is not already attached to an open file
   !
   !> @returns an integer file handle that is not presently associated with any open file
   !---
   function get_unused_file_handle()
         integer :: get_unused_file_handle
         integer :: ierr
         integer, parameter :: minfh = 1001, maxfh = 2001
         logical :: isopen

         do get_unused_file_handle = minfh, maxfh
            inquire(unit=get_unused_file_handle, opened=isopen, iostat=ierr)
            if (ierr /= 0) cycle
            if (.not. isopen) exit
         end do
   end function get_unused_file_handle


   !---
   ! DESCRIPTION:
   !
   !> @brief interpolate n 2d arrays of type real at the same location; return result and derivatives
   !
   !> @details this routine is for multiple arrays because the interpolants only need to be calculated once per set of coordinates
   !> (e.g. [rho,T]), otherwise there is a significant speed hit on the code
   !
   !> @todo even better still may be to just have a routine to calculate the interpolants and then pass them to this routine
   !---
   subroutine generic_interp_2d(xdim, ydim, xco, yco, n, f, myx, myy, ix1, ix2, iy1, iy2, res, dresdx, dresdy)
         integer, intent(in) :: &
               xdim, & !< x-dimension of data
               ydim, & !< y-dimension of data
               ix1, &  !< lower x index
               ix2, &  !< upper x index
               iy1, &  !< lower y index
               iy2, &  !< upper y index
               n       !< number of 2d arrays to interpolate
         real(r8), intent(in) :: &
               xco(xdim), &      !< x-coordinate
               yco(ydim), &      !< y-coordinate
               f(n,xdim,ydim), & !< 2d data of shape (n,xdim,ydim)
               myx, &            !< x-value for which result is desired
               myy               !< y-value for which result is desired
         real(r8) :: x1, x2, y1, y2, dx, dy, xa, xb, ya, yb
         real(r8), dimension(n) :: f11, f21, f12, f22, fA, fB, fC, fD
         real(r8), intent(out), dimension(n) :: &
               res, &    !< result
               dresdx, & !< derivative of result wrt x
               dresdy    !< derivative of result wrt y
         integer :: i

         x1 = xco(ix1); x2 = xco(ix2)
         y1 = yco(iy1); y2 = yco(iy2)
         dx = x2 - x1; dy = y2 - y1

         xb = ( myx - x1 ) / dx; xa = 1._r8 - xb
         yb = ( myy - y1 ) / dy; ya = 1._r8 - yb

         f11(:) = f( :, ix1 , iy1 )
         f21(:) = f( :, ix2 , iy1 )
         f12(:) = f( :, ix1 , iy2 )
         f22(:) = f( :, ix2 , iy2 )

         fA(:)  = f11(:) * ya + f12(:) * yb
         fB(:)  = f12(:) * xa + f22(:) * xb
         fC(:)  = f21(:) * ya + f22(:) * yb
         fD(:)  = f11(:) * xa + f21(:) * xb

         res   (:)   = fA(:) * xa + fC(:) * xb
         dresdx(:)   = ( fC(:) - fA(:) ) / dx
         dresdy(:)   = ( fB(:) - fD(:) ) / dy
   end subroutine generic_interp_2d

   !> sum a 1d array using Kahan summation (compensated summation algorithm)
   function ksum(array)
         real(r8) :: ksum, c, t, y, array(:)
         integer :: i
         c = 0._r8; ksum = 0._r8
         do i = 1, size(array)
            ! subtract compensation
            y = array(i) - c
            ! add contribution
            t = ksum + y
            ! cancel large digits; recover lost part
            c = ( t - ksum ) - y
            ksum = t
         end do
   end function ksum


   !> implements a bisection search for an ascending-ordered r8 type array
   subroutine bsearch_r8( arr, val, idx )
         integer :: &
               idx, &  !< closest index of arr that is less than or equal to val
               siz, &  !< size of array arr
               a,   &  !< lower limit of bisection interval
               b,   &  !< upper limit of bisection interval
               r,   &  !< range ( b - a )
               m       !< mid of bisection interval
         real(r8) :: arr(:), val

         siz = size(arr)

         ! check for out-of-bounds
         if ( val < arr(1) .or. val > arr(siz) ) then
            print *, val, arr(1), arr(siz)
            stop "utils: bsearch_r8: requested value is outside range of array"
         end if

         ! initialise
         a = 1 ; b = siz ; r = b - a ; m = (a + b) / 2

         ! perform binary search
         do while ( val /= arr(m) .and. r > 1 )
            if ( arr(m) > val ) then
               b = m
            else
               a = m
            end if
            r = b - a ; m = (a + b) / 2
         end do

         ! make sure we return the index of arr that is less than val
         if ( arr(m) > val ) m = m - 1

         idx = m
   end subroutine bsearch_r8


   !> This subroutine sets very small negative abundances to a very, very small positive value
   subroutine ypsplus(yps)
         real(r8) :: ypsin(nsp), yps(nsp), yminformal = 1.e-27_r8
         integer :: i

         ypsin(:) = yps(:)
         where ( yps < -yminformal )
            yps = 1.e-99_r8
         end where

         !if (any(ypsin /= yps)) write(*,*) " WARNING: there were large negative abundances in ypsplus"

         where ( yps < yminformal )
            yps = 1.e-99_r8
         end where
   end subroutine ypsplus



end module utils
