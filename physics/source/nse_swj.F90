!> @Author Samuel Jones, Heiko Möller

!> This module solves the NSE constraint equations, formulated as:
!>
!> \f$ \theta = \sum_i X_i (\mu_p,\mu_n) - 1 = 0 \f$
!> 
!> \f$ \phi = \sum_i \frac{Z_i}{A_i} X_i (\mu_p,\mu_n) - Y_\mathrm{e} = 0 \f$
!> 
!> using, optionally, partition functions from the REACLIB database (Cyburt et al 2010) and
!> Coulomb corrections based on Chabrier & Potekhin (1998) and Calder et al (2007). While the
!> equations have been derived long ago, this routine is loosely based on Seitenzahl et al (2009),
!> which is a very nice, and still fairly modern, paper about NSE.
!> 
!> The mass fractions \f$ X_i \f$ of the nuclides in NSE are given by
!> 
!> \f$ X_i = \frac{m_i}{\rho} \frac{g_tot_i}{{\Lambda_i}^3}
!>           \exp\left\{ \frac{Z_i({\mu_p}^\mathrm{kin}+{\mu_p}^\mathrm{coul}) + 
!>           N_i {\mu_n}^\mathrm{kin} - {\mu_i}^\mathrm{coul} + Q_i }{k_B T}  \right\} \f$
!> 
!> where the symbols have their usual meanings, \f$ \Lambda \f$ is the thermal de Broglie
!> wavelength, \f$ g \f$ is the statistical weight \f$ 2J_i + 1 \f$ or the partition function
!> \f$ \sum_i (2J_i + 1) \mathrm{e}^\frac{E_i}{k_B T} \f$ depending on the 
!> \p use_partition_functions switch. \f$ \mu^\mathrm{kin} \f$ and \f$ \mu^\mathrm{coul} \f$ are
!> the kinetic and Coulomb components of the chemical potential and \f$ Q_i \f$ is the nuclear
!> binding energy.
!>
!> The distribution is found for a given \f$ (T,\rho,Y_\mathrm{e}) \f$ using Newton's method
!> to find the proton and neutron chemical potentials:
!> 
!> \f$ \vec{\mu}_{n+1} = \vec{\mu}_n - J^{-1}\vec{V}(\vec{\mu}_n) \f$
!>
!> with
!>
!> \f$ \vec{\mu} \equiv (\mu_p, \mu_n) \f$; \f$ \vec{V} \equiv (\theta, \phi ) \f$
!>
!> and \f$ J \f$ the jacobian matrix \f$ \frac{\partial \vec{V}}{\partial{\vec{\mu}}} \f$
!>
!> The code for the LU decomposition of the Jacobian matrix and backsubstitution step is either
!> from numerical recipes or Frank Timmes or both.
!>

#undef VERBOSE

module nse_swj
   use utils
   use constants
   use array_sizes
   use nuc_data
   use frame_knobs, only: iolevel
   use reaclib, only: exmass, spinpart, anumpart, znumpart, istate, partinter, &
      reaclib_interpolate_partition_functions
   use communication
   use sorting, only: indexx
   implicit none
   private

   !> swtiches
   logical :: &
      use_partition_functions = .true.,  & !< use part. functions from reaclib instead of just GS statistical weights
      use_coulomb_corrections = .true.,  & !< include non-ideal plasma contribution to ion chemical potentials
      write_mu_tables         = .false., & !< rewrite neut and prot kin. chem. pots cache file
      use_mu_tables           = .false., & !< use mu tables for better initial guess
      check_nse_solve         = .false., & !< solve for NSE over a wide range of conditions
      free_nuc_initial        = .false., & !< assume free neutrons and protons based on ye for initial guess
      ni56_initial            = .true.,  & !< assume 100% ni56 for initial guess
      alt_conv                = .true.     !< If true, use theta and phi as a convergence crit, not resid_mu_n and resid_mu_p

   integer, parameter :: &
      CONVERGENCE = 12, &
      BADMATRIX = 13, &
      LNSRCH = 14

   integer, allocatable :: &
      nse2ppn(:), &        !< mapping from nse isotope index to ppn isotope index
      nse2reaclib(:), &    !< mapping from nse isotope index to reaclib isotope index
      iter_cache(:,:,:), & !< iterations cache table
      ierr_cache(:,:,:)    !< error cache table

   integer, public :: &
      ppn2nse(nsp), &     !< mapping from ppn isotope index to nse isotope index
      ppn2reaclib(nsp)     !< mapping from ppn isotope index to reaclib isotope index

   integer :: &
      nse_num_isos, &
      ini56, &
      iprot, &
      ineut, &
      mu_fh, &            !< file handle for chemical potential tables file
      mu_t_dim, &         !< dimensions of chem kin pot (cache) table
      mu_lrho_dim, &      !< dimensions of chem kin pot (cache) table
      mu_ye_dim           !< dimensions of chem kin pot (cache) table

   real(r8) :: &
      nse_jac(2,2), &      !< NSE jacobian matrix
      inv_jac(2,2), &      !< NSE jacobian inverse matrix
      nse_rhs(2), &        !< NSE RHS
      nse_the, &           !< NSE \f$ \theta = \sum_i X_i(\mu_p,\mu_n) - 1 \f$
      nse_phi, &           !< NSE \f$ \phi = \sum_i \frac{Z_i}{A_i} X_i(\mu_p,\mu_n) - Y_\mathrm{e} \f$
      nse_del(2), &        !< \f$ \vec{mu}_{n+1} - \vec{mu}_i \f$
      nse_mu(2), &         !< \f$ [\mu_p^\mathrm{kin}, \mu_n^\mathrm{kin}] \f$
      nse_mu_last(2), &    !< (mu_p, mu_n) from last converged time step
      mu_p_coul, &         !< coulomb correction to proton chemical potential (MeV)
      nnetw(nsp), &        !< neutron numbers of species in ppn network - should go into nuc_data
      anetw(nsp), &        !< mass numbers of species in ppn network
      znetw(nsp), &        !< proton numbers of species in ppn network
      delta_ye, delta_t9, delta_lrho, &
      t9_func,rho_func,ye_func

   ! these are public because they're used by the reverse rate routine
   real(r8), public :: &
      ppn_g(nsp), &        !< statistical weights \f$ g_spin_i = 2J_i + 1 \f$
      ppn_q(nsp), &        !< nuclear binding energies (MeV)
      ppn_m(nsp)           !< nuclear masses (g)

   real(r8), parameter :: &
      mu_rtol               = 1.e-8_r8, & !< N-R rel. tol for convergence on chemical potentials
      thephi_rtol           = 1.e-9_r8, & !< N-R rel. tol for convergence on theta and phi
      y_rtol                = 1.e-5_r8, & !< N-R rel. tol for species
      ye_min                = 0.001_r8, & 
      ye_max                = 0.999_r8, & 
      lrho_min              = 0._r8,    & 
      lrho_max              = 13._r8,   & 
      t9_min                = 1._r8,    & 
      t9_max                = 20._r8

   real(r8), allocatable, dimension(:), public :: &
      g_tot_i, &     !< total statistical weights, with or without accounting for part. function
      mu_i_coul      !< coulomb correction to ion chemical potentials (MeV)

   real(r8), allocatable, dimension(:) :: &
      nse_x, &       !< NSE abundances
      nse_guess, &   !< trial NSE abundances (initial guess for iterative solver)
      g_spin_i, &    !< statistical weights \f$ g_spin_i = 2J_i + 1 \f$
      q_i, &         !< nuclear binding energies (MeV)
      m_i, &         !< nuclear masses (MeV?)
      nse_a, &       !< atomic mass numbers for NSE isotopes
      nse_z, &       !< proton numbers for NSE isotopes
      nse_z_div_a, & !< Z/A for NSE isotopes (cached during init)
      nse_n, &       !< neutron numbers for NSE isotopes
      mu_t9(:), &    !< t9 coordinate of mu (cache) table
      mu_lrho(:), &  !< log10(rho) coordinate of mu (cache) table
      mu_ye(:)       !< ye coordinate of mu (cache) table

   real(r8), allocatable, dimension(:,:,:) :: &
      mu_p_cache(:,:,:), & !< proton kin chem pot cache table
      mu_n_cache(:,:,:)    !< neutron kin chem pot cache table

   character(len=5), allocatable :: nse_zis(:)

   logical :: &
      nse_mask(nsp), &
      first            !< .true. for first time nse solver is called

   character(len=24), parameter ::  mu_fn = '../NPDATA/chemkinpot.dat'

   !> to be destroyed by yoga flame
   common/cnetw/anetw,znetw

   public nse_init, compute_nse, nse_num_isos, coulomb_corrections, partition_functions

contains

   !> the statistical weights and nuclear masses and binding energies should be saved during init
   !> todo nuclear masses should be updated with experimental (e.g. Audi et al), reaclib, FRDM, etc
   !> https://www.oecd-nea.org/dbdata/data/structure.htm#nubase
   !> todo isomers are not treated properly because I don't know where to get the data for them
   !> todo this whole T/F species in networksetup.txt is just so crazy...
   subroutine nse_init()
         implicit none
         integer :: i, j, k, a, z, n, s, ppn_index, idx, ierr
         real(r8) :: yps(nsp), mu_p, mu_n

         if (master) print *, "init nse ... "

         nnetw(:) = anetw(:) - znetw(:)
         ppn_g(:) = ZERO; ppn_m(:) = ZERO; ppn_q(:) = ZERO

         !> statistical weights, masses and binding energies (from mass excesses)
         !> todo for now I take these from reaclib
         do i = 1, nnpartdim
            a = anumpart(i)
            z = znumpart(i)
            s = istate(i)
            s = min(s, 2) ! we only have ground/isom, no thermalised
            n = a - z
            ppn_index = niso( a, z, s )
            if ( ppn_index /= 0 ) then
               ppn_g(ppn_index) = TWO * spinpart(i) + ONE
               ! ^_^ atomic mass in g
               ppn_m(ppn_index) = exmass(i) * mev2erg / clight2 + anetw(ppn_index) * amu
               ! change to nuclear mass in g
               ppn_m(ppn_index) = ppn_m(ppn_index) - z * me_g
               ! ^_^ nuclear binding energy in MeV
               ppn_q(ppn_index) = z * mp_mev + n * mn_mev - ppn_m(ppn_index) * clight2 * erg2mev
               if ( zis(ppn_index) == 'NEUT ' ) ppn_q(ppn_index) = ZERO
               if ( zis(ppn_index) == 'PROT ' ) ppn_q(ppn_index) = ZERO
               ppn2reaclib(ppn_index) = i
            end if
         end do

         !> only solve NSE distribution where we actually want abundances... so stupid
         nse_mask(:) = ( considerisotope .and. anetw /= 0 )

         ! ^_^ isomers other than al*26 we don't have nuclear data for, so ignore for now
         do i = 1, nsp
            if ( nse_mask(i) .and. ( ppn_g(i) == ZERO .or. ppn_m(i) == ZERO ) ) then
               nse_mask(i) = .false.
#ifdef PPN
               print *, ' WARNING: nuclide', zis(i), ' was excluded from NSE state (not enough data)'
#else
               if ( master ) then
                  print *, ' WARNING: nuclide', zis(i), ' was excluded from NSE state (not enough data)'
               endif
#endif
            end if
         end do

         !> set up array sizes and store data
         nse_num_isos = count(nse_mask)
         allocate( nse_x(nse_num_isos), q_i(nse_num_isos), g_spin_i(nse_num_isos), &
               mu_i_coul(nse_num_isos), m_i(nse_num_isos), nse_a(nse_num_isos), &
               nse_z(nse_num_isos), nse_n(nse_num_isos), nse_zis(nse_num_isos), &
               nse2ppn(nse_num_isos), nse2reaclib(nse_num_isos), nse_z_div_a(nse_num_isos), &
               nse_guess(nse_num_isos), g_tot_i(nse_num_isos))

         idx = 0
         do i = 1, nsp
            if ( nse_mask(i) ) then
               idx = idx + 1
               ppn2nse ( i )    = idx
               nse2ppn ( idx )  = i

               q_i         ( idx )  = ppn_q       ( i ) 
               g_spin_i    ( idx )  = ppn_g       ( i ) 
               m_i         ( idx )  = ppn_m       ( i ) 
               nse_a       ( idx )  = anetw       ( i ) 
               nse_z       ( idx )  = znetw       ( i ) 
               nse_n       ( idx )  = nnetw       ( i ) 
               nse_zis     ( idx )  = zis         ( i ) 
               nse2reaclib ( idx )  = ppn2reaclib ( i ) 
            end if
         end do

         nse_z_div_a(:) = nse_z(:) / nse_a(:)

         !> find proton and ni56 coordinate by linear search
         do i = 1, nse_num_isos
            if ( nse_zis(i) == 'NI 56' ) then
               ini56 = i
            end if
            if ( nse_zis(i) == 'PROT ' ) then
               iprot = i
            end if
            if ( nse_zis(i) == 'NEUT ' ) then
               ineut = i
            end if
         end do

         ! caching of neutron and proton chemical potentials for good initial guess
         if ( use_mu_tables ) call load_mu_tables
         if ( write_mu_tables .or. check_nse_solve ) then
            call chemical_potential_preprocessor
         end if

         ! for using last time step's solution as initial gues
         first          = .true.
         nse_mu_last(:) = -9._r8

   end subroutine nse_init


   subroutine compute_nse( t9in, rho, ye, yps, mu_p, mu_n, iter, ierr )
         use errors, only: BAD_YE, BAD_STEP, NSE_CNVRG
         real(r8), intent(in) :: t9in, rho, ye
         real(r8), intent(out) :: yps(nsp), mu_p, mu_n
         integer, intent(out) :: iter, ierr
         real(r8) :: lt9ini, lt9fin, dlt9, lt9,t9
         integer :: i, info, idx, it
         logical :: y_mask(nse_num_isos)
         real(r8) :: coeff
         real(r8) :: t9start = 100.0e0_r8
         integer :: tsteps = 10
         integer, parameter :: max_iter_big = 500
         integer, parameter :: max_iter_small = 100

         ! check input for incosistencies
         if ( rho <= ZERO .or. t9in <= ZERO ) then
            print *, 'problem with input in NSE solver'
            print *, 't9, rho, ye = ', t9, rho, ye
            stop
         end if

         ! bad Ye values can occur during runge-kutta for large time steps. this is
         ! allowed, but an error must be returned so it knows to take a sub-step
         if ( ye <= ZERO .or. ye >= ONE ) then
            ierr = BAD_YE
            return
         end if

         ! try to solve for this t9, rho and ye
         ierr = 999 ; mu_n = ZERO ; mu_p = ZERO ; yps  = ZERO ; it   = 0 ; iter = 0 ; t9   = t9in

         call coulomb_corrections( t9, rho, ye )
         call partition_functions( t9 )
         call initial_guess( mu_p, mu_n, t9, rho, ye )
         nse_guess(:) = nse_distribution( mu_p, mu_n, t9, rho )
         call solve_nse( t9, rho, ye, yps, mu_p, mu_n, it, ierr, max_iter_small )
         iter = iter + it

         select case (ierr)
         case(0)
            ! good solution: set abundances and exit
            nse_x(:) = nse_distribution( mu_p, mu_n, t9, rho )
            yps(:) = 1.e-99_r8
            do i = 1, nse_num_isos
               idx = nse2ppn(i)
               yps(idx) = max( nse_x(i), 1e-99_r8 )
            end do
            ! save solution for next initial guess
            nse_mu_last(:) = nse_mu(:) ; first = .false.
            return
         case(BAD_STEP,NSE_CNVRG)
            ! solving for t9, rho, ye did not work, so continue to temperature descent method
            continue
         case default
            ! something bad happened
            print *, 'problem in NSE solver'
            print *, 't9, rho, ye = ', t9, rho, ye
            stop
         end select

         !> solve indirectly, descending in temperature
         coeff = ( t9start / t9in ) ** ( ONE / (float(tsteps) - ONE) )
         do i = 1, tsteps
            t9 = t9in * coeff ** (tsteps - i)
            call coulomb_corrections( t9, rho, ye )
            call partition_functions( t9 )
            if (i == 1) call initial_guess( mu_p, mu_n, t9, rho, ye )
            nse_guess(:) = nse_distribution( mu_p, mu_n, t9, rho )
            call solve_nse( t9, rho, ye, yps, mu_p, mu_n, it, ierr, max_iter_big )
            iter = iter + it
         enddo

         select case (ierr)
         case(0)
            ! temperature descent method was successful: set abundances and exit
            nse_x(:) = nse_distribution( mu_p, mu_n, t9, rho )
            yps(:) = 1.e-99_r8
            do i = 1, nse_num_isos
               idx = nse2ppn(i)
               yps(idx) = max( nse_x(i), 1e-99_r8 )
            end do
            ! save solution for next initial guess
            nse_mu_last(:) = nse_mu(:) ; first = .false.
            return
         case(NSE_CNVRG)
            ! max. iterations reached in NR/linesearch
            print *, 'nse_swj: newton raphson failed to converge after', max_iter_big, 'iterations'
            print *, 'nse_swj: temperature descent failed'
            print *, 't9, rho, ye = ', t9, rho, ye
            stop
         case(BAD_STEP)
            print *, 'nse_swj: problem in jacobi matrix inversion'
            print *, 't9, rho, ye = ', t9, rho, ye
            stop
         case default
            ! something really bad happened
            print *, 'problem in NSE solver'
            print *, 't9, rho, ye = ', t9, rho, ye
            stop
         end select

   end subroutine compute_nse


   subroutine solve_nse( t9, rho, ye, yps, mu_p, mu_n, it, ierr, max_iter )
         use errors, only: BAD_STEP, NSE_CNVRG
         real(r8), intent(in) :: t9, rho, ye
         real(r8), intent(inout) ::  yps(nsp), mu_p, mu_n
         integer, intent(out) :: it, ierr
         integer, intent(in) :: max_iter
         real(r8) :: perm, resid_p, resid_n, max_resid_y, tti1, tti2
         integer :: i, j, idx, mtx_indx(nse_num_isos), it_lnsrch
         logical :: y_mask(nse_num_isos)
         real(r8) :: gamma
         real(r8) :: nse_the_old,nse_phi_old
         real(r8), dimension(2) :: nse_mu_old
         integer, parameter :: lnsrchmax = 100

         tti1 = wallclocktime() ; nse_mu = [ mu_p, mu_n ] ; nse_x = nse_guess(:) ; it_lnsrch = 0

         !> Newton-Raphson iteration loop
         newton: do i = 1, max_iter
            !> iteration counter
            it = i

            !> store mu old
            nse_mu_old = nse_mu

            !> calculate RHS of NSE constraint equations and Jacobian matrix
            call calculate_rhs( t9, rho, ye ) !based on old nse_x
            !> store old theta and phi for convergence check
            nse_the_old = nse_the !old nse_the
            nse_phi_old = nse_phi !old nse_phi

            !> solve
            call calculate_jac( t9 ) !based on old nse_x
            ierr = 0
            call ludcmp_nse( nse_jac, 2, 2, mtx_indx, perm, ierr )
            if ( ierr /= 0 ) then
               ierr = BAD_STEP
               return
            end if
            call lubksb_nse( nse_jac, 2, 2, mtx_indx, nse_rhs )
            nse_del(:) = nse_rhs(:)

            !> update mu and x and theta and phi
            gamma = ONE
            do j = 1, lnsrchmax
               it_lnsrch = it_lnsrch + 1
               if ( j == lnsrchmax ) then
                  ierr = NSE_CNVRG
                  return
               endif
               nse_mu(:) = nse_mu_old(:) - gamma * nse_del(:)
               nse_x(:) = nse_distribution( nse_mu(1), nse_mu(2), t9, rho )
               ! update nse_the and nse_phi
               call calculate_rhs( t9, rho, ye )
               if (       ( abs(nse_the) > (1.e10_r8 * abs(nse_the_old) + 1.e-10_r8) ) &
                     .or. ( abs(nse_phi) > (1.e10_r8 * abs(nse_phi_old) + 1.e-10_r8) ) ) then
               gamma = gamma * 1.e-2_r8
            else
               exit
            endif
         enddo

#ifdef VERBOSE
         write(*,"(a5, i5,6(es13.5))") 'nse: ', i, mu_p, mu_n, resid_p, resid_n, nse_the, nse_phi, nse_del
#endif

         !> convergence and non-convergence test
         if ( abs(nse_the) <= thephi_rtol .and. abs(nse_phi) <= thephi_rtol ) then
            ! converged
            exit newton
         end if
         if ( i == max_iter ) then
            ! did not converge
            ierr = NSE_CNVRG
            return
         end if

      end do newton

#ifdef VERBOSE
      write(*,*) "newton iterations + linesearch iterations = ", it + it_lnsrch
#endif

      mu_p = nse_mu(1); mu_n = nse_mu(2)
end subroutine solve_nse


!> make an initial guess for \p nse_mu = [ mu_p, mu_n ] either asuming
!> \f$ X(^{56}\mathrm{Fe}) = 1 \f$ and \f$ \mu_p = \mu_n \f$ or using
!> the cache tables
subroutine initial_guess( mu_p, mu_n, t9, rho, ye )
      real(r8), intent(in) :: t9, rho, ye
      real(r8), intent(out) :: mu_p, mu_n
      real(r8) :: lam_ni56, lam_p, lam_n, kT_erg, kT_mev, my_lrho, my_t9, my_ye
      integer :: it9, iye, ilrho

      ! if we have solved previously, use previous time step's solution as initial guess
      if ( .not. first ) then
         nse_mu(:) = nse_mu_last(:)
         mu_p = nse_mu(1) ; mu_n = nse_mu(2)
         return
      end if

      kT_erg = boltzmann * t9 * 1.e9_r8
      kT_mev = kT_erg * erg2mev

      if ( use_mu_tables ) then
         my_t9   = min( max( t9, t9_min ), t9_max )
         my_ye   = min( max( ye, ye_min ), ye_max )
         my_lrho = min( max( log10(rho), lrho_min ), lrho_max )
         it9   = nint( ( my_t9   - t9_min )   / delta_t9 )   + 1
         iye   = nint( ( my_ye   - ye_min )   / delta_ye )   + 1
         ilrho = nint( ( my_lrho - lrho_min ) / delta_lrho ) + 1
         mu_p = mu_p_cache( it9, ilrho, iye )
         mu_n = mu_n_cache( it9, ilrho, iye )
      else if ( free_nuc_initial ) then
         lam_p = ( planck_erg * planck_erg / ( twopi * m_i(iprot) * kT_erg ) ) ** HALF
         lam_n = ( planck_erg * planck_erg / ( twopi * m_i(ineut) * kT_erg ) ) ** HALF
         mu_p  = kT_mev * log( ( ye * rho * lam_p * lam_p * lam_p ) / m_i(iprot) / g_tot_i(iprot) ) &
               &      - q_i(iprot) !+ mu_i_coul(iprot)
         mu_n  = kT_mev * log( ( (ONE - ye) * rho * lam_n * lam_n * lam_n ) / m_i(ineut) / g_tot_i(ineut) ) &
               &      - q_i(ineut) !+ mu_i_coul(ineut)
      else if ( ni56_initial ) then
         lam_ni56 = ( planck_erg * planck_erg / ( twopi * m_i(ini56) * kT_erg ) ) ** HALF
         mu_p = kT_mev * log( ONE * rho * lam_ni56**3 / g_tot_i(ini56) / m_i(ini56) ) - q_i(ini56) + mu_i_coul(ini56) - &
               nse_z(ini56) * mu_p_coul
         mu_p = mu_p / nse_a(ini56)
         mu_n = mu_p
      else
         print *, 'initial guess method not specified, stopping...'
         stop
      end if
#ifdef VERBOSE
      print *, 'initial guess mu_p, mu_n = ', mu_p, mu_n
#endif

end subroutine initial_guess

!> calculate the NSE jacobian matrix
subroutine calculate_jac( t9 )
      real(r8), intent(in) :: t9
      real(r8) :: &
            dthe_dmu_p, &
            dthe_dmu_n, &
            dphi_dmu_p, &
            dphi_dmu_n, &
            dX_dmu_p(nse_num_isos), &
            dX_dmu_n(nse_num_isos), &
            kT_erg, kT_mev
      integer :: i

      kT_erg = boltzmann * t9 * 1.e9_r8
      kT_mev = kT_erg * erg2mev

      dX_dmu_p(:) = nse_x(:) * nse_z(:) / kT_mev
      dX_dmu_n(:) = nse_x(:) * nse_n(:) / kT_mev
      dthe_dmu_p  = sortedsum( dX_dmu_p )
      dthe_dmu_n  = sortedsum( dX_dmu_n )
      dphi_dmu_p  = sortedsum( nse_z_div_a * dX_dmu_p )
      dphi_dmu_n  = sortedsum( nse_z_div_a * dX_dmu_n )

      !> todo check ordering of the jacobian elements
      nse_jac(1,:) = [ dthe_dmu_p, dthe_dmu_n ]
      nse_jac(2,:) = [ dphi_dmu_p, dphi_dmu_n ]
end subroutine calculate_jac

!> calculate the RHS of the NSE equations, \f$ \theta \f$ and \f$ \phi \f$
subroutine calculate_rhs( t9, rho, ye )
      real(r8), intent(in) :: t9, rho, ye
      integer :: i

      nse_the = sortedsum( nse_x )
      nse_phi = sortedsum( nse_z_div_a * nse_x )

      nse_the = nse_the - ONE
      nse_phi = nse_phi - ye

      nse_rhs = [ nse_the, nse_phi ]
end subroutine

!> calculate the NSE distribution from the proton and neutron chemical potentials
function nse_distribution( mu_p, mu_n, t9, rho )
      real(r8), intent(in) :: mu_p         !< prot chemical potential
      real(r8), intent(in) :: mu_n         !< neut chemical potential
      real(r8), intent(in) :: t9
      real(r8), intent(in) :: rho
      real(r8) :: &
            nse_distribution(nse_num_isos), &
            lam(nse_num_isos), &              !< thermal de Broglie wavelengths
            fac(nse_num_isos), &              !< exponential factor
            pre(nse_num_isos), &
            kT_erg, kT_mev
      integer :: i

      kT_erg = boltzmann * t9 * 1.e9_r8
      kT_mev = kT_erg * erg2mev
      lam(:) = ( planck_erg * planck_erg / ( twopi * m_i(:) * kT_erg ) ) ** HALF      
      fac(:) = nse_z(:) * ( mu_p + mu_p_coul ) + nse_n(:) * mu_n - mu_i_coul(:) + q_i(:)
      pre(:) = m_i(:) / rho * g_tot_i(:) / ( lam(:) * lam(:) * lam(:) )
      nse_distribution(:) = pre(:) * exp( fac(:) / kT_mev )
end function nse_distribution


subroutine coulomb_corrections( t9, rho, ye )
      real(r8), intent(in) :: t9, rho, ye
      real(r8) :: kT_mev, kT_erg, &
            Gam_i(nse_num_isos), & !< ion coulomb coupling parameter
            Gam_e, &               !< electron coulomb coupling parameter
            ae                     !< electron sphere radius
      real(r8), parameter :: &
            a1 = -0.9052_r8, &                                   !< coulomb corrections constant
            a2 =  0.6322_r8, &                                   !< coulomb corrections constant
            a3 = - HALF * (THREE ** HALF) - a1 / (a2 ** HALF)    !< coulomb corrections constant
      integer :: i

      kT_erg = boltzmann * t9 * 1.e9_r8
      kT_mev = kT_erg * erg2mev
      if ( use_coulomb_corrections ) then
         ae = ( FOUR * pi * rho * ye * avogadro / THREE ) ** (-THIRD)
         Gam_e = qe * qe / kT_erg / ae
         do i = 1, nse_num_isos
            Gam_i(i) = Gam_e * nse_z(i) ** FIVE_THIRDS
            mu_i_coul(i) = a1 * ( &
                  ( Gam_i(i)*(a2 + Gam_i(i)) )**HALF - &
                  a2 * log( (Gam_i(i)/a2)**HALF + (ONE + (Gam_i(i)/a2))**HALF ) &
                  ) + &
                  TWO * a3 * ( Gam_i(i)**HALF - atan( Gam_i(i)**HALF ) )
            mu_i_coul(i) = mu_i_coul(i) * kT_mev
         end do
         mu_p_coul = mu_i_coul(iprot)
      else
         mu_p_coul    = ZERO
         mu_i_coul(:) = ZERO
      end if
end subroutine coulomb_corrections


!> get partition functions from reaclib and save for nse
subroutine partition_functions( t9 )
      real(r8), intent(in) :: t9
      integer :: i

      if ( use_partition_functions ) then
         call reaclib_interpolate_partition_functions( t9 )         
         do i = 1, nse_num_isos
            g_tot_i(i) = partinter( nse2reaclib(i) ) * g_spin_i(i)
         end do
      else
         g_tot_i(:) = g_spin_i(:)
      end if
end subroutine partition_functions



subroutine chemical_potential_preprocessor
      integer :: i, j, k, error_count, iters, ierr
      real(r8) :: yps(nsp), t9, rho, ye, mu_n, mu_p

      if ( write_mu_tables ) then
         write(*,*) 'Writing chemical potential tables'
         if ( use_mu_tables ) stop 'nse_swj: cannot use and write mu tables simultaneously'
      end if

      mu_t_dim    = 41 !< dimensions of chem kin pot (cache) table
      mu_lrho_dim = 14 !< dimensions of chem kin pot (cache) table
      mu_ye_dim   = 21 !< dimensions of chem kin pot (cache) table

      call allocate_mu_tables


#ifndef PPN
      if ( master ) then
#endif

         mu_fh = get_unused_file_handle()
         if ( write_mu_tables ) open( mu_fh, file = mu_fn, form = "unformatted", access = "stream" )
         delta_t9 = ( t9_max - t9_min ) / ( mu_t_dim - 1 )
         delta_ye = ( ye_max - ye_min ) / ( mu_ye_dim - 1 )
         delta_lrho = ( lrho_max - lrho_min ) / ( mu_lrho_dim - 1 )
         mu_t9(:)   = (/ ( t9_min   + (i - 1) * delta_t9   , i = 1, mu_t_dim    ) /)
         mu_lrho(:) = (/ ( lrho_min + (i - 1) * delta_lrho , i = 1, mu_lrho_dim ) /)
         mu_ye(:)   = (/ ( ye_min   + (i - 1) * delta_ye   , i = 1, mu_ye_dim   ) /)

         error_count = 0 ; iter_cache(:,:,:)  = 0

         do i = 1, mu_t_dim
            do j = 1, mu_lrho_dim
               do k = 1, mu_ye_dim
                  t9 = mu_t9(i); rho = 10._r8**mu_lrho(j); ye = mu_ye(k)
                  first = .true.
                  if ( write_mu_tables ) then
                     call compute_nse( t9, rho, ye, yps, &
                           mu_p_cache(i,j,k), mu_n_cache(i,j,k), & 
                           iter_cache(i,j,k), ierr_cache(i,j,k) )
                  else 
                     call compute_nse( t9, rho, ye, yps, mu_p, mu_n, &
                           iter_cache(i,j,k), ierr_cache(i,j,k) )
                  end if
                  if ( ierr_cache(i,j,k) /= 0 ) error_count = error_count + 1
                  print *, t9, rho, ye, iter_cache(i,j,k), error_count,sum(yps), mu_p_cache(i,j,k), mu_n_cache(i,j,k)
               end do 
            end do
         end do

         if ( write_mu_tables ) then
            write( mu_fh ) mu_t_dim, mu_lrho_dim, mu_ye_dim
            write( mu_fh ) mu_t9(:), mu_ye(:), mu_lrho(:)
            write( mu_fh ) delta_t9, delta_lrho, delta_ye
            write( mu_fh ) mu_p_cache(1:mu_t_dim,1:mu_lrho_dim,1:mu_ye_dim)
            write( mu_fh ) mu_n_cache(1:mu_t_dim,1:mu_lrho_dim,1:mu_ye_dim)
            write( mu_fh ) iter_cache(1:mu_t_dim,1:mu_lrho_dim,1:mu_ye_dim)
            write( mu_fh ) ierr_cache(1:mu_t_dim,1:mu_lrho_dim,1:mu_ye_dim)
            close( mu_fh )
            stop 'mu tables for NSE written to NPDATA'
         end if
         if ( check_nse_solve ) stop 'nse routine tested.'
#ifndef PPN
      end if
#endif
end subroutine chemical_potential_preprocessor


subroutine load_mu_tables
      if ( write_mu_tables ) stop 'nse_swj: cannot write and use mu tables simultaneously'
#ifndef PPN
      if ( master ) then
#endif
         call check_file_exists( mu_fn )
         mu_fh = get_unused_file_handle()
         open( mu_fh, file = mu_fn, form = "unformatted", access = "stream" )
         read( mu_fh ) mu_t_dim, mu_lrho_dim, mu_ye_dim
         call allocate_mu_tables
         read( mu_fh ) mu_t9(:), mu_ye(:), mu_lrho(:)
         read( mu_fh ) delta_t9, delta_lrho, delta_ye
         read( mu_fh ) mu_p_cache(:,:,:)
         read( mu_fh ) mu_n_cache(:,:,:)
         close( mu_fh )
#ifndef PPN
      end if
      call nse_broadcasts
#endif
end subroutine load_mu_tables


! allocate memory for chemical potential tables
subroutine allocate_mu_tables
      allocate ( mu_t9(mu_t_dim), mu_lrho(mu_lrho_dim), mu_ye(mu_ye_dim), &
            mu_n_cache(mu_t_dim,mu_lrho_dim,mu_ye_dim), &
            mu_p_cache(mu_t_dim,mu_lrho_dim,mu_ye_dim), &
            iter_cache(mu_t_dim,mu_lrho_dim,mu_ye_dim), &
            ierr_cache(mu_t_dim,mu_lrho_dim,mu_ye_dim) )
      iter_cache = 0    ; ierr_cache = 0
      mu_p_cache = ZERO ; mu_n_cache = ZERO
      mu_t9      = ZERO ; mu_lrho    = ZERO ; mu_ye = ZERO
end subroutine allocate_mu_tables



#ifndef PPN
subroutine nse_broadcasts
      call broadcast(mu_t_dim) ; call broadcast(mu_lrho_dim) ; call broadcast(mu_ye_dim)
      if ( .not. allocated(mu_t9) ) call allocate_mu_tables
      call broadcast(mu_t9)
      call broadcast(mu_ye)
      call broadcast(mu_lrho)
      call broadcast(delta_t9)
      call broadcast(delta_lrho)
      call broadcast(delta_ye)
      call broadcast(mu_p_cache)
      call broadcast(mu_n_cache)
end subroutine nse_broadcasts
#endif




subroutine ludcmp_nse(a,n,np,indx,d,ierr)
      implicit none
      save

      !..given the matrix a(n,n), with physical dimensions a(np,ap) this routine
      !..replaces a by the lu decompostion of a row-wise permutation of itself.
      !..input are a,n,np. output is a, indx which records the row
      !..permutations effected by the partial pivoting, and d which is 1 if
      !..the number of interchanges is even, -1 if odd.
      !..use routine lubksb to solve a system of linear equations.
      !..
      !..nmax is the largest expected value of n
      ! ^_^ ierr = 1 for singular matrix; ierr=0 and AOK

      !..declare
      integer          n,np,indx(np),nmax,i,j,k,imax,ierr
      parameter        (nmax=500)
      double precision a(np,np),d,tiny,vv(nmax),aamax,sum,dum
      parameter        (tiny=1.0d-20)


      !..vv stores the implicit scaling of each row
      !..loop over the rows to get the scaling information
      d = 1.0d0
      do i=1,n
         aamax = 0.0d0
         do j=1,n
            if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
         enddo
         if (aamax  ==  ZERO) then
            ierr = 1
            return
            !stop 'singular matrix in ludcmp'
         end if
         vv(i) = 1.0d0/aamax
      enddo

      !..for each column apply crouts method; see equation 2.3.12
      do j=1,n
         do i=1,j-1
            sum = a(i,j)
            do k=1,i-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
         enddo

         !..find the largest pivot element
         aamax = 0.0d0
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum .ge. aamax) then
               imax  = i
               aamax = dum
            end if
         enddo

         !..if we need to interchange rows
         if (j .ne. imax) then
            do k=1,n
               dum       = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k)    = dum
            enddo
            d          = -d
            vv(imax)   = vv(j)
         end if

         !..divide by the pivot element
         indx(j) = imax
         if (a(j,j)  ==  ZERO) a(j,j) = tiny
         if (j .ne. n) then
            dum = 1.0d0/a(j,j)
            do i=j+1,n
               a(i,j) = a(i,j)*dum
            enddo
         end if

         !..and go back for another column of crouts method
      enddo
end subroutine ludcmp_nse


subroutine lubksb_nse(a,n,np,indx,b)
      implicit none
      save

      !..solves a set of n linear equations ax=b. a is input in its lu decomposition
      !..form, determined by the routine above ludcmp. indx is input as the
      !..permutation vector also returned by ludcmp. b is input as the right hand
      !..side vector and returns with the solution vector x.
      !..a,n ans np are not modified by this routine and thus can be left in place
      !..for successive calls (i.e matrix inversion)

      !..declare
      integer           n,np,indx(np),i,ii,j,ll
      double precision  a(np,np),b(np),sum

      !..when ii is > 0, ii becomes the index of the first nonzero element of b
      !..this is forward substitution of equation 2.3.6, and unscamble in place
      ii = 0
      do i=1,n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii .ne. 0) then
            do j=ii,i-1
               sum = sum - a(i,j) * b(j)
            enddo

            !..nonzero element was found, so dos the sums in the loop above
         else if (sum .ne. ZERO) then
            ii  = i
         end if
         b(i) = sum
      enddo

      !..back substitution equation 2.3.7
      do i = n,1,-1
         sum = b(i)
         if (i .lt. n) then
            do j=i+1,n
               sum = sum - a(i,j) * b(j)
            enddo
         end if
         b(i) = sum/a(i,i)
      enddo
end subroutine lubksb_nse

!    subroutine matinv(a,n,y)
!       implicit none
! 
!       real(r8) :: a(n,n),y(n,n),d
!       integer ::  i,j,n,indx(n),ierr
!       y=0.
!       do i=1,n
!         y(i,i)=1.
!       end do
!       call ludcmp_nse(a,n,n,indx,d,ierr)
!       do j=1,n
!         call lubksb_nse(a,n,n,indx,y(1,j))
!       end do
! 
!    end subroutine matinv

end module nse_swj
