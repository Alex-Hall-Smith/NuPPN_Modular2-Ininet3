! TODO replace num_tracers with ntracers or something consistent with data names
module traj_data
   use utils
   use constants
   use trajectories
   use abundances
   use frame_knobs, only: ini_filename, restart, interpolation_type
   use nuc_data, only: epsi
   use communication
   implicit none
   private

   integer(i4), public :: nsteps_max, nsteps_min, num_tracers, num_done, stride, pnum
   integer(i4), private :: fh, ifirstcycle, ilastcycle, nxt_pnum, pxt_ptcl, abn
   integer(i4), allocatable, public :: &
         tracers_done(:), &  !< for restarting
         pnums(:)

   ! io
   character(len=256), public :: datdir
   character(len=256) :: filein
   character(len=256), parameter, public :: outfile = './output/abu.dat'
   character(len=256), parameter, public :: donefile = './output/num_done.dat'
   character(len=256), parameter, public :: failfile = './output/failures.txt'

   ! file handles
   integer(i4), public :: TRAJ_INPUT, ABU_OUTPUT, FAIL_OUTPUT
   integer(i4), private :: RESTART_LIST, RESTART_NUM

   integer(i4), public :: nspecies
   integer(i8), public :: last_write_offset, last_read_offset

   ! trajectory object
   type ( trajectory ), public :: t

   real(r8), private :: trajmenclosed

   public :: traj_data_init, read_trajectory, dump_trajectory

   integer  :: one_int
 
contains


   subroutine traj_data_init()
         implicit none

         if (master) print *, " init tracer particle data ... "
         ! file handle assignment
         TRAJ_INPUT   = get_unused_file_handle()
         ABU_OUTPUT   = TRAJ_INPUT   + 1
         FAIL_OUTPUT  = ABU_OUTPUT   + 1

         call traj_data_bin_init()
   end subroutine traj_data_init


   subroutine open_bin_input()
         integer(i4) :: ierr

         ! open the file
         open(newunit = TRAJ_INPUT, file = trim(filein), form = "unformatted", access = "stream", &
               action = "read", status = "old")

   end subroutine open_bin_input


   subroutine close_bin_input()
         close(TRAJ_INPUT)
   end subroutine close_bin_input


   subroutine traj_data_bin_init()
         use utils, only: bigend
         use frame_knobs, only: tppnp_input_filename
         integer :: i

         filein = trim(datdir) // tppnp_input_filename

         ! set offset for input file
         last_read_offset = 0
         
         call check_file_exists( filein )
         call open_bin_input()

         ! read the total number of trajectories and min/max number of time steps from the header
         read( TRAJ_INPUT ) num_tracers
         read( TRAJ_INPUT ) nsteps_min
         read( TRAJ_INPUT ) nsteps_max

         ! set needle
         inquire(TRAJ_INPUT,pos=last_read_offset)
         call close_bin_input()

   end subroutine traj_data_bin_init


   !> read one trajectory from the input file
   subroutine read_trajectory()
         implicit none

         call read_one_trajectory

         if ( restart ) then
            ! keep reading trajectories until we find one that has not already been computed
            do while ( any( tracers_done(:) == pnum ) )
               call read_one_trajectory()
            end do
         end if

         t%temp(:) = t%temp(:) / 1.e9_r8

         ! set up interpolant for trajectory
         ! NB - this needs to be done by the slave upon receipt, so currently this won't work in
         ! tppnp.
         !call create_traj_interpolants(t, interpolation_type)
   end subroutine read_trajectory

   subroutine read_one_trajectory
         call read_bin_trajectory
   end subroutine read_one_trajectory


   !> reads a trajectory and its initial abundances
   !> from the lsb/msb fortran binary trajectory file
   subroutine read_bin_trajectory()
         integer(i4) :: i, offset
         integer(i4) :: entry_size !< number of bytes in a particle data set

         ! deallocate time evolution arrays
         if ( allocated(t%time) ) then
            deallocate(t%time) ; deallocate(t%temp) ; deallocate(t%dens)
         end if

         ! deallocate initial abundances arrays
         if ( allocated( t%xini_x ) ) then
            deallocate( t%xini_x ) ; deallocate( t%xini_a ) ; deallocate( t%xini_z )
            deallocate( t%xini_zis )
         end if

         ! open the input file at the correct position for the next tracer particle
         call open_bin_input()

         ! read time-evolution data, starting from where we left off
         read( TRAJ_INPUT, pos = last_read_offset ) pnum, t%ntimesteps
         allocate( t%time(t%ntimesteps), t%temp(t%ntimesteps), t%dens(t%ntimesteps) )

         ! read mass of tracer particle and initial and final positions
         read( TRAJ_INPUT ) t%mass
         read( TRAJ_INPUT ) t%posini(:)
         read( TRAJ_INPUT ) t%posfin(:)

         read( TRAJ_INPUT ) t%time(:)
         read( TRAJ_INPUT ) t%temp(:)
         read( TRAJ_INPUT ) t%dens(:)

         ! read initial abundance data
         read( TRAJ_INPUT ) pnum, t%xini_num_isos
         allocate( t%xini_x(t%xini_num_isos), t%xini_zis(t%xini_num_isos) , &
               t%xini_a(t%xini_num_isos), t%xini_z(t%xini_num_isos) )
         read( TRAJ_INPUT ) t%xini_a
         read( TRAJ_INPUT ) t%xini_z
         read( TRAJ_INPUT ) t%xini_zis
         read( TRAJ_INPUT ) t%xini_x

         ! set variables from the abundances module for abundance initialisation
         num_iso_in = t % xini_num_isos
         iso_in ( 1 : num_iso_in ) = t % xini_zis(:)
         x_in   ( 1 : num_iso_in ) = max( t % xini_x(:), 1.e-99_r8 )

         inquire(TRAJ_INPUT,pos=last_read_offset)

         call close_bin_input()
   end subroutine read_bin_trajectory 

   !> dump a trajectory in ascii to the file "dump.trajectory" for human inspection
   subroutine dump_trajectory(t)
      implicit none
      type(trajectory) :: t
      integer(i4) :: i, fh

      open(newunit = fh, file = "dump.trajectory", action = "write")
      write(fh,*) t%ntimesteps
      do i = 1, t%ntimesteps
         write(fh,"(3(es23.15))") t%time(i), t%temp(i), t%dens(i)
      end do

      write(fh,*) t%xini_num_isos
      do i = 1, t%xini_num_isos
         write(fh,"(i5,i5,2x,a5,es23.15)") t%xini_a(i), t%xini_z(i), t%xini_zis(i), t%xini_x(i)
      end do

      stop "trajectory dumped to dump.trajectory"
   end subroutine dump_trajectory

end module traj_data

