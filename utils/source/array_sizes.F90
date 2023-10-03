!> @brief provides general array sizes needed across the code
!> @details the array sizes are mostly defined during compilation using preprocessor variable definitions based on entries in the
!isotopedatabase.txt file and the ppn_frame.input and ppn_physics.input files
!
!> @todo for now nre1 is given from beta decay tables; in case bigger tables are introduced, nre1 must be updated!
!> @todo the array sizes should ideally be compilation-independent, but we wanted to eradicate multiple parameter files

module array_sizes
   implicit none
   public

   integer, parameter :: partdim   = 24   !< number of partition functions values given vs temperature
   integer, parameter :: nre1      = 2500 !< max number of rates in netgen tables
   integer, parameter :: iAtdim    = 85   !< max atomic number in reaclib network (index_reaclib = 1,2)
   integer, parameter :: i282dim   = 282  !< max mass number in reaclib network (index_reaclib = 1,2)
   integer, parameter :: iCfdim    = 98   !< max atomic number in reaclib network (index_reaclib = 3)
   integer, parameter :: i325dim   = 325  !< max mass number in reaclib network (index_reaclib = 3)
   integer, parameter :: ngrid     = 100  !< size of netgen tables

#ifdef pNNN
   integer, parameter :: nsp = pNNN + pNVCP + 20 !< max number of species in network
   integer, parameter :: nre = nsp * 15          !< max number of reactions in network
#if defined(PPN) || defined(TPPNP)
   integer, parameter :: msl = 1            !< number of mass shells (1 in PPN)
#else
   integer, parameter :: msl     = pmsl     !< number of mass shells (1 in PPN)
   integer, parameter :: nrefmax = pnrefmax !< number of sepcies to refine grid on (AMR depreciated)
   integer, parameter :: gfdim   = pgfdim
#endif
   
#if pIDX_RCLB == 0
   integer, parameter :: ndim_rl   = 61721 !< number of rates in reaclib
   integer, parameter :: nnpartdim = 5427  !< number of isotopes with partition functions in reaclib
#elif pIDX_RCLB == 1
   integer, parameter :: ndim_rl   = 74313 !< number of rates in reaclib
!   integer, parameter :: nnpartdim = 5492  !< number of isotopes with partition functions in reaclib
   integer, parameter :: nnpartdim = 5456  !< number of isotopes with partition functions in reaclib
! RJS -- INDEX_REACLIB=1 reads a different partition function data file than =2, which has fewer lines
#elif pIDX_RCLB == 2
   integer, parameter :: ndim_rl   = 74313 !< number of rates in reaclib
!   integer, parameter :: nnpartdim = 5492  !< number of isotopes with partition functions in reaclib
   integer, parameter :: nnpartdim = 5497   ! number of isotopes with partition functions in NEW reaclib test 27/04/22
#elif pIDX_RCLB == 3
   integer, parameter :: ndim_rl   = 77953 !< number of rates in reaclib
   integer, parameter :: nnpartdim = 6789  !< number of isotopes with partition functions in reaclib
#else
#error "INDEX_REACLIB not defined"
#endif
   ! (pIDX_RCLB)
#else
   include 'parameter.inc'
#endif
   !(pNNN)
end module array_sizes

