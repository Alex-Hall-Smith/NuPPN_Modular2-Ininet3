!> This module contains diagnostic properties of the various solvers.
module solver_diagnostics
   use utils
   implicit none
   private

   integer, public :: &
         iter, &            !< iteration number
         nsubt, &           !< number of sub-timesteps
         titers, &          !< total iterations of Newton-Raphson method for entire time step integration
         nvar1

   real(r8), public :: &
         time_matrix_inv, & !< time taken for matrix inversion
         ntime, &           !< total time spent in one Newton-Raphson iteration
         ntime_last, &      !< time taken for one Newton-Raphson iteration
         time_jacobian

end module solver_diagnostics
