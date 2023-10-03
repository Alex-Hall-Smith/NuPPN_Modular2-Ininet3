module trajectories
   use utils
   use interp
   use constants
   use akima, only: create_akima_interpolant, akima_interpolant, interp_akima
   implicit none
   private

   character(len=256), public :: trajectory_fn
   integer :: traj_fh

   type, public :: trajectory
      real(r8), allocatable :: &
            time(:), & !< time / s
            temp(:), & !< temperature / K
            dens(:), & !< density / g/cc
            xini_x(:)  !< initial composition (mass fractions)
      real(r8) :: &
            mass, &       !< mass of tracer particle
            tpeak, &      !< peak temperature
            rhopeak, &    !< peak density
            posini(3), &  !< (initial) position: x,y,z or just enclosed mass
            posfin(3)     !< (final) position: x,y,z or just enclosed mass
      integer(i4), allocatable  :: &
            xini_a(:), &   !< atomic mass number of isotopes for which initial composition is given
            xini_z(:)      !< proton number of isotopes for which initial composition is given
      integer(i4) :: &
            ntimesteps, &    !< number of time steps
            xini_num_isos, & !< number of isotopes for which initial composition is given
            tdim, &          !< dimension of simulation from which simulation was taken (1/2/3)
            pnum
      character(len=5), allocatable :: &
            xini_zis(:)      !< names of sotopes for which initial composition is given
      character(len=3) :: &
            temp_unit, & !< 'KKK', 'T8K' or 'T9K'
            dens_unit, & !< 'CGS' or 'LOG', presumably meaning g/cc or log10(g/cc)
            time_unit    !< 'SEC', 'YRS', 'DTY' or 'DTS'
      logical :: timestep_given
      integer(i4) :: interpolation_type ! linear (1) or akima (3)
      type(interpolant) :: li_temp, li_dens ! linear interpolants in temperature and density
      type(akima_interpolant) :: ai_temp, ai_dens ! akima interpolant in temperature and densitt
   end type trajectory

   public ppn_trajectory_init, read_ppn_trajectory, create_traj_interpolants, interpolate_trajectory

contains

   subroutine ppn_trajectory_init(t, interpolation_type)
         type(trajectory) :: t
         integer(i4) :: interpolation_type

         call read_ppn_trajectory(t)

         ! ^_^ ensure trajectories are in seconds, GK and g/cc

         select case(t%time_unit)
         case('YRS','DTY')
            t%time(:) = t%time(:) * yrs2sec
         end select

         select case(t%dens_unit)
         case('LOG')
            t%dens(:) = 10._r8 ** t%dens(:)
         end select

         select case(t%temp_unit)
         case('T8K')
            t%temp(:) = t%temp / 10._r8
         case('KKK')
            t%temp(:) = t%temp / 1.e9_r8
         end select

         select case(t%time_unit)
         case('DTY','DTS')
            t%timestep_given = .true.
         case default
            t%timestep_given = .false.
         end select

         call create_traj_interpolants(t, interpolation_type)

   end subroutine ppn_trajectory_init


   subroutine read_ppn_trajectory(t)
         integer :: i, ierr
         type(trajectory) :: t
         traj_fh = get_unused_file_handle()

         ! ^_^ nasty way to get number of timesteps
         open(traj_fh, file = trajectory_fn)
         do i = 1, 7
            read( traj_fh, * )
         end do
         do i = 0, 99999999
            read( traj_fh, *, iostat = ierr )
            if ( ierr < 0 ) then
               t%ntimesteps = i
               exit
            end if
         end do

         allocate( t%temp(t%ntimesteps) )
         allocate( t%time(t%ntimesteps) )
         allocate( t%dens(t%ntimesteps) )

         ! ^_^ read trajectory
         rewind(traj_fh)
         do i = 1, 3
            read( traj_fh, * )
         end do
         read( traj_fh, '(10x,A3)' ) t%time_unit
         read( traj_fh, '(10x,A3)' ) t%temp_unit
         read( traj_fh, '(10x,A3)' ) t%dens_unit
         read( traj_fh, * )
         do i = 1, t%ntimesteps
            read( traj_fh, * ) t%time(i), t%temp(i), t%dens(i)
         end do
         close(traj_fh)
   end subroutine read_ppn_trajectory


   subroutine create_traj_interpolants(t, interpolation_type)
         type(trajectory) :: t
         integer(i4) :: interpolation_type
         if (interpolation_type == 0) return ! no interpolant needed
         if (t%timestep_given) stop "trajectory interpolation currently not possible in timestep given mode"
         t%interpolation_type = interpolation_type
         select case (t%interpolation_type)
         case(1) ! linear
            call create_1d_interpolant( t%ntimesteps, t%time, t%temp, t%li_temp )
            call create_1d_interpolant( t%ntimesteps, t%dens, t%temp, t%li_dens )
         case(3) ! akima
            call create_akima_interpolant( t%time, t%temp, t%ai_temp )
            call create_akima_interpolant( t%time, t%dens, t%ai_dens )
         end select
   end subroutine create_traj_interpolants


   subroutine interpolate_trajectory(t, time, temp, dens)
         type(trajectory) :: t
         real(r8), intent(in) :: time
         real(r8), intent(out) :: temp, dens
         select case(t%interpolation_type)
         case(1)
            call interpolate_1d( t%li_temp, time, temp )
            call interpolate_1d( t%li_dens, time, dens )
         case(3)
            call interp_akima( t%ai_temp, time, temp )
            call interp_akima( t%ai_dens, time, dens )
         end select
   end subroutine interpolate_trajectory
end module trajectories


