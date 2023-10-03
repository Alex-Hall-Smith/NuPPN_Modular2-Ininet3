! simple module to provide an interface to the integrator for decaying composition for some time
module decays
   use utils
   use array_sizes
   use physics_knobs, only: decay_time
   use nuc_data, only: considerisotope
   implicit none

   real(r8) :: anetw(nsp), znetw(nsp)
   common/cnetw/anetw,znetw

   public :: do_decay

   contains

      subroutine do_decay(nvar, yps, nvrel)
         use neutrinos
         use solver, only: integrate_network
         implicit none
         integer(i4) :: nvar, nvrel, ierr
         real(r8) :: yps(nsp), ye
         real(r8), parameter :: t9 = 0.01_r8, rho = 10._r8
         type(neutrino) :: nu

         call init_one_neutrino(nu)
         call calculate_ye(yps, anetw, znetw, ye, considerisotope)

         ierr = 0
         call integrate_network( &
            nvar, yps, t9, t9, rho, rho, ye, decay_time, nvrel, nu, ierr)

         if (ierr /= 0) stop "decay was unsuccessful"
      end subroutine do_decay
end module decays
