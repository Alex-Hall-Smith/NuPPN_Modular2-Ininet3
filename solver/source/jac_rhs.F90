module jac_rhs
   use array_sizes
   use utils
   use constants
   use solver_knobs
   use nuc_data
   use sparsekit_module, only: getelm, coocsr, csrcsc
   use solver_diagnostics
   use frame_knobs, only: iolevel
   use communication
   implicit none
   private

   real(r8) :: &
         fact(0:4)       !< factorials

   type(mtx_coo), public :: jac_coo         !< Coordinate-format sparse jacobian matrix for largest possible network
   type(mtx_coo), public :: jac_cpy         !< Matrix copy

   type(mtx_csr), public :: jac_csr         !< Row-compressed (CSR) format sparse jacobian matrix for largest possible network

   type(mtx_csr), public :: reduced_jac_csr !< Row-compressed (CSR) jacobian matrix for current iteration of dynamic network
   type(mtx_csc), public :: reduced_jac_csc !< Col-compressed (CSC) jacobian matrix for current iteration of dynamic network
   type(mtx_coo), public :: reduced_jac_coo !< Coordinate-format jacobian matrix for current iteration of dynamic network

   real(r8), allocatable, public :: reduced_rhs(:), reduced_y_new(:), reduced_y_cur(:), reduced_dxdt(:)

   real(r8), public :: &
         dxdt(nsp), & !< time derivative of mass fractions
         d2xdtdt(nsp) !< temperature derivative of dxdt

   integer :: &
         dns2csr_map(8,nre), & !< for each reaction, the 8 indices of its possible contributions to the sparse jacobian matrix
         in1in1, in2in1, out1in1, out2in1, in1in2, in2in2, out1in2, out2in2, &
         iamgamma, &     !< index of gamma in network
         nnz, &          !< number of non-zeros in full jacobian matrix
         reduced_nnz

   public :: calculate_dxdt, calculate_jacobian, makemap, jac_rhs_init, nvar1, dns2csr, &
         free_reduced_jacobian, free_reduced_rhs, reduce_jacobian, reduce_rhs, reduce_dxdt, &
         reduce_rhs_bd, ddecompac, nnz, reduced_nnz, convert_reduced_jac_to_csc
#ifndef PPN
   public :: jac_rhs_broadcasts
#endif

contains

   subroutine jac_rhs_init()
         integer :: i

         ! find gamma
         do i = 1, nsp
            if ( zis(i) == 'OOOOO' ) then
               iamgamma = i
               exit
            end if
         end do

         ! factorials
         fact(0) = ONE
         fact(1) = ONE
         fact(2) = HALF
         fact(3) = SIXTH
         fact(4) = ONE / 24._r8

         ! ^_^ mapping between dense and CSR jacobian matrices,
         !     to ensure contiguous memory access in index mapping array:
         in1in1  =           1
         in2in1  =  in1in1 + 1
         out1in1 =  in2in1 + 1
         out2in1 = out1in1 + 1
         in1in2  = out2in1 + 1
         in2in2  =  in1in2 + 1
         out1in2 =  in2in2 + 1
         out2in2 = out1in2 + 1

         nnz = 0
   end subroutine jac_rhs_init


   subroutine free_reduced_jacobian
         implicit none
         call free_csr_matrix( reduced_jac_csr )
         call free_coo_matrix( reduced_jac_coo )
   end subroutine free_reduced_jacobian


   subroutine free_reduced_rhs
         implicit none
         deallocate( reduced_rhs )
         if ( allocated( reduced_y_cur ) ) deallocate( reduced_y_cur )
         if ( allocated( reduced_y_new ) ) deallocate( reduced_y_new )
   end subroutine free_reduced_rhs


   !> @brief Reduce RHS and species arrays for current problem size
   subroutine reduce_rhs( nvar1, mask, y_new, y_cur, rhs )
         implicit none
         integer :: nvar1
         real(r8) :: y_new(nsp), y_cur(nsp), rhs(nsp)
         logical :: mask(nsp)

         allocate( reduced_rhs(nvar1), reduced_y_new(nvar1), reduced_y_cur(nvar1) )
         reduced_rhs   = pack ( rhs    , mask ) 
         reduced_y_new = pack ( y_new  , mask ) 
         reduced_y_cur = pack ( y_cur  , mask ) 
   end subroutine reduce_rhs


   subroutine reduce_rhs_bd( nvar1, mask, rhs )
         implicit none
         integer(i4) :: nvar1
         real(r8), intent(in) :: rhs(nsp)
         logical :: mask(nsp)

         allocate( reduced_rhs(nvar1) )
         reduced_rhs   = pack ( rhs, mask )
   end subroutine reduce_rhs_bd


   subroutine reduce_dxdt( nvar1, mask )
         implicit none
         integer :: nvar1
         logical :: mask(nsp)

         allocate( reduced_dxdt(nvar1) )
         reduced_dxdt = pack( dxdt, mask )
   end subroutine reduce_dxdt


   !> @brief Add identity matrix (1 in diagonal) to jacobian matrix, then
   !> reduce to problem size (remove columns/rows with only a 1 in diagonal) then convert
   !> from coordinate format to CSR format
   subroutine reduce_jacobian( mask, nvar )
         implicit none
         integer :: i, ice, nvar, n_ice, iexpand(nsp)
         logical :: mask(nsp)

         !> flag remaining zeros in jacobian (not sure why they are there actually)
         do i = 1, size(jac_coo%val)
            if ( ( ( abs(jac_coo%val(i)) == ZERO ) .or. jac_coo%val(i) == -ONE ) .and. jac_coo%row(i) /= -1) then
               if ( iolevel >= 2 ) print *, 'full: zero at ', i, ':', jac_coo%val(i), jac_coo%row(i)
               jac_coo%row(i) = -1
               if ( iolevel >= 2) print *, '      row ammended: ', jac_coo%row(i)
            end if
         end do

         ! save codec for dynamic network
         nvar1 = 0
         iexpand(:) = -1
         do i = 1, nsp
            if ( .not. mask(i) ) cycle
            nvar1 = nvar1 + 1
            iexpand(i) = nvar1

            ! ^_^ add 1 to diagonal only for elements we are including in the dynamic-size solve
            ice = dns2csr(i,i)
            jac_coo % row(ice) = i
            jac_coo % col(ice) = i
            jac_coo % val(ice) = jac_coo % val(ice) + ONE
         end do

         n_ice = count ( jac_coo % row /= -1 )

         ! ^_^ pack coordinate-format jacobian matrix (squeeze out zeros)
         call allocate_coo_matrix( reduced_jac_coo, n_ice )
         reduced_jac_coo % row(:) = pack( jac_coo % row , jac_coo % row /= -1 )
         reduced_jac_coo % col(:) = pack( jac_coo % col , jac_coo % row /= -1 )
         reduced_jac_coo % val(:) = pack( jac_coo % val , jac_coo % row /= -1 )

         ! ^_^ correct coordinate indices from full matrix (1-nsp) to reduced matrix (1-nvar1)
         do i = 1, n_ice
            reduced_jac_coo % row(i) = iexpand(reduced_jac_coo % row(i))
            reduced_jac_coo % col(i) = iexpand(reduced_jac_coo % col(i))
         end do

         if ( any( abs( reduced_jac_coo%val ) == ZERO ) ) then
            do i = 1, n_ice
               if ( abs ( reduced_jac_coo%val(i) ) == ZERO ) print *, 'reduced: zero at ', i
            end do
            stop 'zero in sparse matrix'
         end if

         ! ^_^ convert to CSR format
         call allocate_csr_matrix( reduced_jac_csr, n_ice, nvar1 )
         call coocsr( nvar1, n_ice, &
               reduced_jac_coo % val, &
               reduced_jac_coo % row, &
               reduced_jac_coo % col, &
               reduced_jac_csr % val, &
               reduced_jac_csr % col, &
               reduced_jac_csr % row )

         ! record number of non-zeros in reduced jacobian matrix
         reduced_nnz = n_ice

   end subroutine reduce_jacobian


   subroutine convert_reduced_jac_to_csc()
         implicit none
         integer :: nnz, n
         nnz = size(reduced_jac_csr%val)
         n   = size(reduced_jac_csr%row) - 1
         allocate( reduced_jac_csc%val (   nnz ) )
         allocate( reduced_jac_csc%col ( n + 1 ) )
         allocate( reduced_jac_csc%row (   nnz ) )

         call csrcsc( n, 1, 1, &
               reduced_jac_csr%val, reduced_jac_csr%col, reduced_jac_csr%row, &
               reduced_jac_csc%val, reduced_jac_csc%row, reduced_jac_csc%col )

   end subroutine


   !> @brief Calculates the time derivatives of the mass fractions \f$ \frac{\partial X}{\partial t} \f$
   subroutine calculate_dxdt( x_in, dt_est )
         use array_sizes, only: nsp, nre
         use nuc_data, only: considerisotope, considerreaction
         use rates, only: k1, k2, k3, k4, k5, k6, k7, k8, v, dvdt
         use frame_knobs, only : iplot_flux_option
         use physics_knobs, only: getderivs
         implicit none
         real(r8), intent(in) :: x_in(nsp) !< isotopic mass fractions
         real(r8), intent(out) :: dt_est  !< time step estimate based on fastest relevant time scale
         real(r8) :: anetw(nsp), znetw(nsp), xn(nsp), vx, flux_to_print(nre), &
               dx1dt, dx3dt, dx5dt, dx7dt, timescales(4)
         integer :: i
         common /cnetw/ anetw, znetw
         common /flux/ flux_to_print

         dxdt(:) = ZERO; xn(:) = ZERO ; dt_est = 1e99_r8 

         ! ^_^ calculate molar abundances
         where (considerisotope .and. anetw /= 0)
            xn(:) = x_in(:) / anetw(:)
         end where

         ! set molar abundance of gamma to 1
         xn(iamgamma) = ONE

         do i = 1, nre
            if ( .not. considerreaction(i) ) cycle

            vx = fact(k2(i)) * v(i) * xn(k1(i)) ** k2(i) * xn(k3(i)) ! (production/destruction flux.)
            if (iplot_flux_option == 1) then
               flux_to_print(i) = vx
            end if

            dx1dt = - k2(i) * vx ; dx3dt = -vx ; dx5dt = k6(i) * vx ; dx7dt = k8(i) * vx

            dxdt(k1(i)) = dxdt(k1(i)) + dx1dt
            dxdt(k3(i)) = dxdt(k3(i)) + dx3dt
            dxdt(k5(i)) = dxdt(k5(i)) + dx5dt
            dxdt(k7(i)) = dxdt(k7(i)) + dx7dt

            ! time step estimator based on abundance and flux
            !timescales = abs( [ &
                  !max(xn(k1(i)),cyminformal) / max(dx1dt,SMALL), &
                  !max(xn(k3(i)),cyminformal) / max(dx3dt,SMALL), &
                  !max(xn(k5(i)),cyminformal) / max(dx5dt,SMALL), &
                  !max(xn(k7(i)),cyminformal) / max(dx7dt,SMALL) &
                  !] )

            !dt_est = min( dt_est, minval( timescales, mask = (timescales > ZERO) ) )

         end do

         ! convert to mass fractions
         where (considerisotope)
            dxdt(:) = dxdt(:) * anetw(:)
         end where

         ! set molar abundance of gamma and dxdt(igamma) back to 0
         xn(iamgamma) = ZERO
         dxdt(iamgamma) = ZERO

   end subroutine calculate_dxdt


   !> This subroutine calculates the Jacobian matrix of the nuclear network
   subroutine calculate_jacobian( dzeit, mask, x_in )
         use nuc_data
         use rates
         use frame_knobs
         implicit none

         real(r8), intent(in) :: x_in(nsp)
         real(r8), dimension(nsp) :: anetw, znetw, xn
         real(r8) :: dzeit, dfdN, stem, val, dfdNmin, &
               dfin1dy1, dfin2dy1, dfout1dy1, dfout2dy1, dfin1dy2, dfin2dy2, dfout1dy2, dfout2dy2
         integer(i4) :: i, l
         logical :: mask(nsp)
         common /cnetw/ anetw, znetw                                

         select case( irdn )
         case(1)
            stop "irdn1 is no longer supported"
         case(2)
            continue
         case default
            stop "jacobn: choice of irdn not found"
         end select

         associate(jac => jac_coo % val, row => jac_coo % row, col => jac_coo % col)

            mask(:) = .false.
            dfdNmin = cyminformal / dzeit
            jac(:) = 1.e-99_r8 ; row(:) = -1 ; col(:) = -1

            ! calculate molar abundances
            where ( considerisotope(:) .and. anetw(:) > ZERO )
               xn(:) = x_in(:) / anetw(:)
            end where

            ! set molar abundance of gamma to 1
            xn(iamgamma) = ONE

            do l = 1, nre 
               if ( .not. considerreaction(l) ) cycle
               if ( abs(v(l)) / dzeit <= min_rate_ov_dt ) cycle

               ! ^_^ in1, in2 are the reactants, of which there are num_in1 and num_in2
               !     out1, out2 are the products, of which there are num_out1, num_out2
               associate( &
                        in1   => k1(l)        , in2   => k3(l)        , &
                        out1  => k5(l)        , out2  => k7(l)        , &
                        Nin1  => k2(l)        , Nin2  => k4(l)        , &
                        Nout1 => k6(l)        , Nout2 => k8(l)        , &
                        ain1  => anetw(k1(l)) , ain2  => anetw(k3(l)) , &
                        aout1 => anetw(k5(l)) , aout2 => anetw(k7(l)) , &
                        y1    => xn(k1(l))    , y2    => xn(k3(l)), &
                        ! > jacobian coordinates
                        ice1  => dns2csr_map(  in1in1 , l ), &
                        ice2  => dns2csr_map(  in2in1 , l ), &
                        ice3  => dns2csr_map(  in1in2 , l ), &
                        ice4  => dns2csr_map( out1in1 , l ), &
                        ice5  => dns2csr_map( out2in1 , l ), &
                        ice6  => dns2csr_map(  in2in2 , l ), &
                        ice7  => dns2csr_map( out1in2 , l ), &
                        ice8  => dns2csr_map( out2in2 , l ) &
                        )

                  ! ^_^ ( f == dydt )
                  stem = fact(Nin1) * fact(Nin2) * v(l)
                  dfdN = stem * y1 ** Nin1 * y2

                  ! ^_^ select species to consider in network solve
                  ! (i.e. size of jacobian) based on dydt*deltat
                  if ( abs(dfdN) > dfdNmin ) then
                     mask(  in1 )  = .true.;  mask(  in2 )  = .true.
                     mask( out1 )  = .true.;  mask( out2 )  = .true.
                  end if

                  stem = Nin1 * fact(Nin1) * fact(Nin2) * y2 ** Nin2 * y1 ** (Nin1 - 1) * v(l)
                  dfin1dy1  = - Nin1  * stem; dfin2dy1  = - Nin2  * stem
                  dfout1dy1 =   Nout1 * stem; dfout2dy1 =   Nout2 * stem

                  stem = Nin2 * fact(Nin1) * fact(Nin2) * y1 ** Nin1 * y2 ** (Nin2 - 1) * v(l)
                  dfin1dy2  = - Nin1  * stem; dfin2dy2  = - Nin2  * stem
                  dfout1dy2 =   Nout1 * stem; dfout2dy2 =   Nout2 * stem

                  ! store jacobian in COO format

                  if ( in1 /= iamgamma ) then
                     jac(ice1) = jac(ice1) + dfin1dy1 * ain1 / ain1
                  end if

                  if ( in1 /= iamgamma .and. in2 /= iamgamma ) then
                     jac(ice2) = jac(ice2) + dfin2dy1 * ain2 / ain1
                     jac(ice3) = jac(ice3) + dfin1dy2 * ain1 / ain2
                  end if

                  if ( in1 /= iamgamma .and. out1 /= iamgamma ) then
                     jac(ice4) = jac(ice4) + dfout1dy1 * aout1 / ain1
                  end if

                  if ( in1 /= iamgamma .and. out2 /= iamgamma ) then
                     jac(ice5) = jac(ice5) + dfout2dy1 * aout2 / ain1
                  end if

                  if ( in2 /= iamgamma ) then
                     jac(ice6) = jac(ice6) + dfin2dy2 * ain2 / ain2
                  end if

                  if ( in2 /= iamgamma .and. out1 /= iamgamma ) then
                     jac(ice7) = jac(ice7) + dfout1dy2 * aout1 / ain2
                  end if

                  if ( in2 /= iamgamma .and. out2 /= iamgamma ) then
                     jac(ice8) = jac(ice8) + dfout2dy2 * aout2 / ain2
                  end if

                  ! here p, n and alpha are excluded from the network solve if they are only
                  ! important for capture reactions, i.e. (a,g), (g,p) (n,g), ... etc,
                  ! though the impact of these reactions on other species (i.e. the main
                  ! reactant and product of such reactions) are still important

                  if ( in2  /= iamgamma .and. .not. mask(in2)  ) cycle
                  if ( out1 /= iamgamma .and. .not. mask(out1) ) cycle

                  ! the matrix is now written in coordinate format, where the i and j indices
                  ! of each non-zero element of the jacobian is stored in the arrays xx and col_coo.
                  ! Although all reactions have contributed to the elements of the jacobian, we
                  ! only perform the network solve for the part of the matrix that matters
                  ! (according to the cyminformal criteria of the dynamic network).
                  if ( mask(in1) .and. mask(out2) ) then
                     if ( in1  /=  iamgamma ) then
                        row(ice1) = in1; col(ice1) = in1
                     end if

                     if (in1  /=  iamgamma .and. in2  /=  iamgamma) then
                        row(ice2) = in2; col(ice2) = in1 
                        row(ice3) = in1; col(ice3) = in2
                     end if

                     if ( in1  /=  iamgamma .and. out1  /=  iamgamma ) then
                        row(ice4) = out1; col(ice4) = in1
                     end if

                     if ( in1  /=  iamgamma .and. out2  /=  iamgamma ) then
                        row(ice5) = out2; col(ice5) = in1
                     end if

                     if ( in2  /=  iamgamma ) then
                        row(ice6) = in2; col(ice6) = in2
                     end if

                     if ( in2  /=  iamgamma .and. out1  /=  iamgamma) then
                        row(ice7) = out1; col(ice7) = in2
                     end if

                     if ( in2  /=  iamgamma .and. out2  /=  iamgamma) then
                        row(ice8) = out2; col(ice8) = in2
                     end if
                  end if

               end associate
            end do

            ! set molar abundance of gamma back to 0 and switch off mask
            xn(iamgamma) = ZERO; mask(iamgamma) = .false.

            do l = 1, size(row)
               if (jac(l) == ZERO) then
                  if (row(l) /= col(l)) then
                     print *, 'zero at', row(l), col(l), zis(row(l)), zis(col(l))
                     stop 'non-diagonal zero in sparse matrix'
                  end if
               end if
            end do

         end associate

   end subroutine calculate_jacobian

   !---
   ! DESCRIPTION:
   !> @brief build the map between the jacobian matrix non-zeros (i.e. where we have reactions) and the sparse matrix format, and
   !> allocate the largest possible sparse matrix in coordinate and row-compressed formats
   subroutine makemap()
         use nuc_data
         use rates, only: k1, k3, k5, k7

         integer :: i, j, cnt, in1, in2, out1, out2
         integer, allocatable :: mtx2sps_2d_temp(:,:)

         allocate( mtx2sps_2d_temp( nsp,nsp ) )
         mtx2sps_2d_temp(:,:) = -1
         iamgamma = ispe('OOOOO')

         ! set potential non-zeros to 1
         do i = 1, nre
            if ( .not. considerreaction(i) ) cycle
            in1     = k1(i) ;  in2     = k3(i) ; out1     = k5(i) ; out2     = k7(i)
            if (  in1  ==  0 ) stop 'k1 is zero in makemap'
            if (  in2  ==  0 ) stop 'k3 is zero in makemap'
            if ( out1  ==  0 ) stop 'k5 is zero in makemap'
            if ( out2  ==  0 ) stop 'k7 is zero in makemap'

            mtx2sps_2d_temp (  in1 , in1 ) = 1
            mtx2sps_2d_temp (  in2 , in1 ) = 1
            mtx2sps_2d_temp ( out1 , in1 ) = 1
            mtx2sps_2d_temp ( out2 , in1 ) = 1
            mtx2sps_2d_temp (  in1 , in2 ) = 1
            mtx2sps_2d_temp (  in2 , in2 ) = 1
            mtx2sps_2d_temp ( out1 , in2 ) = 1
            mtx2sps_2d_temp ( out2 , in2 ) = 1
         end do

         ! Ensure 1 in the diagonal:
         do i = 1, nsp
            if ( considerisotope(i) ) mtx2sps_2d_temp( i, i ) = 1
         end do

         nnz = count( mtx2sps_2d_temp == 1 )
         ! allocate full jacobian matrix in coordinate format
         call allocate_coo_matrix( jac_coo, nnz )
         call allocate_coo_matrix( jac_cpy, nnz )

         ! tag them in row-order, numbering them sequentially
         cnt = 0
         do i = 1, nsp
            do j = 1, nsp
               if ( mtx2sps_2d_temp(i,j) /= -1 ) then
                  cnt = cnt + 1
                  mtx2sps_2d_temp(i,j) = cnt
                  jac_coo % row(cnt) = i ; jac_coo % col(cnt) = j ; jac_coo % val(cnt) = cnt
               end if
            end do
         end do
         deallocate( mtx2sps_2d_temp )

         ! ^_^ allocate and write full jacobian matrix in CSR format
         call allocate_csr_matrix( jac_csr, nnz, nsp )
         call coocsr( nsp, nnz, &
               jac_coo % val, &
               jac_coo % row, &
               jac_coo % col, &
               jac_csr % val, &
               jac_csr % col, &
               jac_csr % row &
               )
            if (master) write(*,*) 'set up sparse matrix'

         ! ^_^ set up mapping. entries in the sparse jacobian should be indexed in this order
         !     to ensure contiguous memory access in index mapping array:
         do i = 1, nre
            if ( .not. considerreaction(i) ) cycle
            in1     = k1(i) ;  in2     = k3(i) ; out1     = k5(i) ; out2     = k7(i)
            dns2csr_map ( in1in1 , i ) = dns2csr( in1 , in1 )
            dns2csr_map ( in2in1 , i ) = dns2csr( in2 , in1 )
            dns2csr_map ( out1in1, i ) = dns2csr( out1, in1 )
            dns2csr_map ( out2in1, i ) = dns2csr( out2, in1 )
            dns2csr_map ( in1in2 , i ) = dns2csr( in1 , in2 )
            dns2csr_map ( in2in2 , i ) = dns2csr( in2 , in2 )
            dns2csr_map ( out1in2, i ) = dns2csr( out1, in2 )
            dns2csr_map ( out2in2, i ) = dns2csr( out2, in2 )
         end do
   end subroutine makemap


   !---
   ! DESCRIPTION
   !
   !> @brief Translate dense Jacobian matrix coordinates into CSR-format Jacobian matrix coordinate
   !> @returns index in full CSR jacobian matrix
   !---
   function dns2csr( row, col )
         integer, intent(in) :: &
               row, & !< row in full dense jacobian matrix
               col    !< column in full dense jacobian matrix
         integer :: dns2csr, row_ptr
         real(r8) :: val

         val = getelm( row, col, &
               jac_csr % val, &
               jac_csr % col, &
               jac_csr % row, &
               dns2csr, .true. )
   end function dns2csr


#ifndef PPN
   subroutine jac_rhs_broadcasts
         call broadcast(nnz)
         call broadcast(dns2csr_map)
         if ( .not. master ) then
            call allocate_csr_matrix( jac_csr, nnz, nsp )
            call allocate_coo_matrix( jac_coo, nnz )
            call allocate_coo_matrix( jac_cpy, nnz )
         end if
         call broadcast(jac_csr % val)
         call broadcast(jac_csr % col)
         call broadcast(jac_csr % row)
         call broadcast(jac_coo % val)
         call broadcast(jac_coo % col)
         call broadcast(jac_coo % row)
         call broadcast(iamgamma)
   end subroutine jac_rhs_broadcasts
#endif


   subroutine ddecompac( mask, dcy, ddy, m, nsp, nvar, dy )

      integer :: m, nsp, nvar, icons, i, j, lc, L
      real(r8) :: ddy(m*nsp), dy(m*nsp), dcy(m*nvar)
      logical :: mask(nsp)

      icons = 0
      do i = 1, nsp
         if (mask(i)) then
            icons = icons + 1
            do j = 1, M 
               LC = (j-1) * nvar
               L  = (j-1) * nsp
               dy(L+i) = dcy(LC+icons)
            end do
         else
            do j = 1, M 
               L = (j-1) * nsp
               dy(L+i) = ddy(L+i) 
            end do
         end if
      end do

      if (icons /= nvar) stop "ddecompac: error in decompac"
   end subroutine ddecompac


end module jac_rhs
