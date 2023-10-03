!> @author S. W. Jones
!
!> @brief Provides supprt and interfaces for MPI communication protocols
!
!> @details The MPPNP and TPPNP frames are parallelised using MPI (Message Passing Interface). This module provides several
!> interfaces and subroutines, utilities, etc to facilitate using MPI protocol in the parallel application

module communication
   use utils, only: netgen_data, r8, i4
   implicit none
   private

#ifdef PPN

   logical, parameter, public :: master = .true.
   integer, parameter, public :: ipid = 0

   !> Generic dump interface
   interface dump
      module procedure dump_r8_array, dump_i4_array, dump_l_array, dump_ch_array
   end interface

   public dump

contains

   !> dump an logical array to a file with the processor's name
   subroutine dump_l_array(array, halt)
         logical :: array(:), halt
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,2x,L)") i, array(i)
         end do
         
         close(fh)
         if (halt) stop
   end subroutine dump_l_array


   !> dump an i4 array to a file with the processor's name
   subroutine dump_i4_array(array, halt)
         integer(i4) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,i15)") i, array(i)
         end do
         
         close(fh)
         if (halt) stop
   end subroutine dump_i4_array


   !> dump an r8 array to a file with the processor's name
   subroutine dump_r8_array(array, halt)
         real(r8) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,es25.15e3)") i, array(i)
         end do
         
         close(fh)
         if (halt) stop
   end subroutine dump_r8_array

   !> dump character array to a file with the processor's name
   subroutine dump_ch_array(array, halt)
         character(len=*) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,2x,a)") i, array(i)
         end do
         
         close(fh)
         if (halt) stop
   end subroutine dump_ch_array


#else

   include "mpif.h"
   integer, public :: num_procs                !< number of tasks running
   integer, public :: ipid                     !< process/task id
   integer, public :: imaster                  !< master task id (=0)
   logical, public :: master                   !< true if process is master
   logical, public :: slave                    !< true if process is slave
   integer, allocatable, public :: disp_mpi(:)
   integer, allocatable, public :: rcount(:)
   integer :: ierr                             !< error flag

   !> Generic broadcast interface
   interface broadcast
      module procedure broadcast_si_0d, broadcast_si_1d, broadcast_si_2d, broadcast_si_3d, broadcast_si_4d, &
            broadcast_sr_0d, broadcast_sr_1d, broadcast_sr_2d, broadcast_sr_3d, &
            broadcast_dr_0d, broadcast_dr_1d, broadcast_dr_2d, broadcast_dr_3d, &
            broadcast_l_0d,  broadcast_l_1d
   end interface

   !> Generic interface to MPI_allreduce with MPI_MAX
   interface allreduce_max
      module procedure allreduce_max_dr_0d, allreduce_max_dr_2d, allreduce_max_si_0d
   end interface

   !> Generic interface to MPI_reduce  with MPI_MAX
   interface reduce_max
      module procedure reduce_max_dr_1d, reduce_max_si_1d
   end interface

   !> Generic interface to MPI_allreduce with MPI_SUM
   interface allreduce_sum
      module procedure allreduce_sum_dr_2d, allreduce_sum_dr_1d, allreduce_sum_si_1d
   end interface

   !> Generic dump interface
   interface dump
      module procedure dump_r8_array, dump_i4_array, dump_l_array, dump_ch_array
   end interface


   public mpi_comm_world, mpi_tag, mpi_source, mpi_any_tag, mpi_integer, mpi_real, mpi_integer4, &
         mpi_real8, mpi_double_precision, mpi_status_size, mpi_any_source, mpi_logical, &
         mpi_packed, comm_init, mpi_integer8, broadcast_netgen_data, allreduce_max, allreduce_sum, &
         reduce_max, terminate_process, shutdown, broadcast, broadcast_ch_0d, broadcast_ch_arr, &
         dump

contains

   !> MPI communication initialization routine
   subroutine comm_init
         call mpi_init(ierr)
         call mpi_comm_rank(mpi_comm_world, ipid, ierr)
         call mpi_comm_size(mpi_comm_world, num_procs, ierr)
         imaster = 0
         master  = ( ipid == imaster )
         slave   = ( ipid /= imaster )
         allocate( disp_mpi(num_procs), rcount(num_procs) )
   end subroutine comm_init

   subroutine terminate_process(taskid)
         integer(i4) :: taskid
         call mpi_send(1, 0, mpi_integer, taskid, 0, mpi_comm_world, ierr)
   end subroutine terminate_process

   !> dump an logical array to a file with the processor's name
   subroutine dump_l_array(array,halt)
         logical :: array(:), halt
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,2x,L)") i, array(i)
         end do
         
         close(fh)

         if (halt) then
            call mpi_barrier(mpi_comm_world)
            stop
         end if
   end subroutine dump_l_array


   !> dump an i4 array to a file with the processor's name
   subroutine dump_i4_array(array,halt)
         integer(i4) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,i15)") i, array(i)
         end do
         
         close(fh)

         if (halt) then
            call mpi_barrier(mpi_comm_world)
            stop
         end if
   end subroutine dump_i4_array


   !> dump an r8 array to a file with the processor's name
   subroutine dump_r8_array(array,halt)
         real(r8) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,es23.15)") i, array(i)
         end do
         
         close(fh)

         if (halt) then
            call mpi_barrier(mpi_comm_world)
            stop
         end if
   end subroutine dump_r8_array

   !> dump character array to a file with the processor's name
   subroutine dump_ch_array(array,halt)
         character(len=*) :: array(:)
         integer :: i, fh
         character(len=10) :: fn
         character(len=5) :: num, prefix = "dump_"
         logical :: halt

         write(num,"(i5)") ipid
         fn = prefix // trim(adjustl(num))

         open(newunit = fh, file = fn)
         do i = 1, size(array)
            write(fh,"(i7,2x,a)") i, array(i)
         end do
         
         close(fh)

         if (halt) then
            call mpi_barrier(mpi_comm_world)
            stop
         end if
   end subroutine dump_ch_array



   !> double real reduce max (1d)
   subroutine reduce_max_dr_1d(val)
         double precision :: val(:)
         integer :: ierr
         if (master) then
            call mpi_reduce( mpi_in_place, val, size(val), mpi_double_precision, mpi_max, imaster, mpi_comm_world, ierr )
         else
            call mpi_reduce( val, val, size(val), mpi_double_precision, mpi_max, imaster, mpi_comm_world, ierr )
         end if
   end subroutine reduce_max_dr_1d

   !> single integer reduce max (1d)
   subroutine reduce_max_si_1d(val)
         integer :: val(:), ierr
         if (master) then
            call mpi_reduce( mpi_in_place, val, size(val), mpi_integer, mpi_max, imaster, mpi_comm_world, ierr )
         else
            call mpi_reduce( val, val, size(val), mpi_integer, mpi_max, imaster, mpi_comm_world, ierr )
         end if
   end subroutine reduce_max_si_1d

   !> single integer allreduce max
   subroutine allreduce_max_si_0d(val)
         integer :: val, ierr
         call mpi_allreduce( mpi_in_place, val, 1, mpi_integer, mpi_max, mpi_comm_world, ierr )
   end subroutine allreduce_max_si_0d

   !> double real allreduce max
   subroutine allreduce_max_dr_0d(val)
         double precision :: val
         integer :: ierr
         call mpi_allreduce( mpi_in_place, val, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr )
   end subroutine allreduce_max_dr_0d

   !> double real allreduce max (2d)
   subroutine allreduce_max_dr_2d(val)
         double precision :: val(:,:)
         integer :: ierr
         call mpi_allreduce( mpi_in_place, val, size(val), mpi_double_precision, mpi_max, mpi_comm_world, ierr )
   end subroutine allreduce_max_dr_2d

   !> double real allreduce sum (2d)
   subroutine allreduce_sum_dr_2d(val)
         double precision :: val(:,:)
         integer :: ierr
         call mpi_allreduce( mpi_in_place, val, size(val), mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
   end subroutine allreduce_sum_dr_2d

   !> double real allreduce sum (1d)
   subroutine allreduce_sum_dr_1d(val)
         double precision :: val(:)
         integer :: ierr
         call mpi_allreduce( mpi_in_place, val, size(val), mpi_double_precision, mpi_sum, mpi_comm_world, ierr )
   end subroutine allreduce_sum_dr_1d

   !> single integer allreduce sum (1d)
   subroutine allreduce_sum_si_1d(val)
         integer :: val(:)
         integer :: ierr
         call mpi_allreduce( mpi_in_place, val, size(val), mpi_integer, mpi_sum, mpi_comm_world, ierr )
   end subroutine allreduce_sum_si_1d

   !> 0d character broadcasting
   subroutine broadcast_ch_0d(val)
         character(len=*) :: val
         integer :: length
         length = len(val)
         call mpi_bcast( val, length, mpi_character, 0, mpi_comm_world, ierr )
   end subroutine broadcast_ch_0d

   !> variable character array broadcasting
   subroutine broadcast_ch_arr(val)
         character(len=*) :: val(:)
         integer :: length
         length = len(val(1))
         call mpi_bcast( val, length*size(val), mpi_character, 0, mpi_comm_world, ierr )
   end subroutine broadcast_ch_arr

   !> single precision integer
   subroutine broadcast_si_0d(val)
         integer :: val
         call mpi_bcast( val, 1, mpi_integer, 0, mpi_comm_world, ierr )
   end subroutine broadcast_si_0d

   !> single integer array broadcasting (1d)
   subroutine broadcast_si_1d(val)
         integer :: val(:)
         call mpi_bcast( val, size(val), mpi_integer, 0, mpi_comm_world, ierr )
   end subroutine broadcast_si_1d

   !> single integer array broadcasting (2d)
   subroutine broadcast_si_2d(val)
         integer :: val(:,:)
         call mpi_bcast( val, size(val), mpi_integer, 0, mpi_comm_world, ierr )
   end subroutine broadcast_si_2d

   !> single integer array broadcasting (3d)
   subroutine broadcast_si_3d(val)
         integer :: val(:,:,:)
         call mpi_bcast( val, size(val), mpi_integer, 0, mpi_comm_world, ierr )
   end subroutine broadcast_si_3d

   !> single integer array broadcasting (4d)
   subroutine broadcast_si_4d(val)
         integer :: val(:,:,:,:)
         call mpi_bcast( val, size(val), mpi_integer, 0, mpi_comm_world, ierr )
   end subroutine broadcast_si_4d

   !> single precision real
   subroutine broadcast_sr_0d(val)
         real :: val
         call mpi_bcast( val, 1, mpi_real, 0, mpi_comm_world, ierr )
   end subroutine broadcast_sr_0d

   !> single real array broadcasting (1d)
   subroutine broadcast_sr_1d(val)
         real :: val(:)
         call mpi_bcast( val, size(val), mpi_real, 0, mpi_comm_world, ierr )
   end subroutine broadcast_sr_1d

   !> single real array broadcasting (2d)
   subroutine broadcast_sr_2d(val)
         real :: val(:,:)
         call mpi_bcast( val, size(val), mpi_real, 0, mpi_comm_world, ierr )
   end subroutine broadcast_sr_2d

   !> single real array broadcasting (3d)
   subroutine broadcast_sr_3d(val)
         real :: val(:,:,:)
         call mpi_bcast( val, size(val), mpi_real, 0, mpi_comm_world, ierr )
   end subroutine broadcast_sr_3d

   !> double precision real
   subroutine broadcast_dr_0d(val)
         double precision :: val
         call mpi_bcast( val, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
   end subroutine broadcast_dr_0d

   !> double precision real array (1d)
   subroutine broadcast_dr_1d(val)
         double precision :: val(:)
         call mpi_bcast( val, size(val), mpi_double_precision, 0, mpi_comm_world, ierr )
   end subroutine broadcast_dr_1d

   !> double precision real array (2d)
   subroutine broadcast_dr_2d(val)
         double precision :: val(:,:)
         call mpi_bcast( val, size(val), mpi_double_precision, 0, mpi_comm_world, ierr )
   end subroutine broadcast_dr_2d

   !> double precision real array (3d)
   subroutine broadcast_dr_3d(val)
         double precision :: val(:,:,:)
         call mpi_bcast( val, size(val), mpi_double_precision, 0, mpi_comm_world, ierr )
   end subroutine broadcast_dr_3d

   !> logical broadcast
   subroutine broadcast_l_0d(val)
         logical :: val
         call mpi_bcast( val, 1, mpi_logical, 0, mpi_comm_world, ierr )
   end subroutine broadcast_l_0d

   !> logical array broadcast (1d)
   subroutine broadcast_l_1d(val)
         logical :: val(:)
         call mpi_bcast( val, size(val), mpi_logical, 0, mpi_comm_world, ierr )
   end subroutine broadcast_l_1d

   !> broadcast derived NETGEN data type
   subroutine broadcast_netgen_data(val)
         type(netgen_data) :: val

         call broadcast(val% tgrid)
         call broadcast(val% l10t8grid)
         call broadcast(val% aNegrid)
         call broadcast(val% vgrid)
         call broadcast(val% l10vgrid)
         call broadcast(val% flag)
         call broadcast(val% Qrad)
         call broadcast(val% Qnu)
         call broadcast(val% kgrid)
         call broadcast(val% ireac)
         call broadcast(val% a1)
         call broadcast(val% a2)
         call broadcast(val% a3)
         call broadcast(val% a4)
         call broadcast(val% n1)
         call broadcast(val% n2)
         call broadcast(val% n3)
         call broadcast(val% n4)
         call broadcast(val% ndata)
         call broadcast(val% zt)
         call broadcast(val% zf)
         call broadcast(val% at)
         call broadcast(val% af)
         call broadcast_ch_arr(val% z1)
         call broadcast_ch_arr(val% z2)
         call broadcast_ch_arr(val% z3)
         call broadcast_ch_arr(val% z4)
         call broadcast_ch_arr(val% elt)
         call broadcast_ch_arr(val% elf)

   end subroutine broadcast_netgen_data


   subroutine shutdown(msg)
         character(len=*) msg
         integer :: ierr
         if ( master ) print *, 'shutdown: ', msg
         call mpi_finalize( ierr )
         stop
   end subroutine shutdown
#endif
!PPN

end module communication















