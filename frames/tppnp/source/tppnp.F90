program tppnp

   !     88P'888'Y88 888 88e  888 88e  Y88b Y88 888 88e  
   !     P'  888  'Y 888 888D 888 888D  Y88b Y8 888 888D 
   !         888     888 88^  888 88^  b Y88b Y 888 88^  
   !         888     888      888      8b Y88b  888      
   !         888     888      888      88b Y88b 888      
   !
   !
   !     TPPNP: this is the trajectory post-processing network/parallel 
   !     frame, which is part of the NuPPN code suite developed and maintained
   !     by the NuGrid collaboration (www.nugridstars.org)
   !     (c) Sam Jones (2016-)

   use array_sizes
   use traj_data
   use frame_knobs
   use solver_knobs
   use physics_knobs
   use solver, only: solver_init, solver_broadcasts, integrate_network
   use backward_euler
   use nse_solver
   use solver_diagnostics
   use jac_rhs
   use nuc_data
   use rates, only: rates_init, k1, k2, k3, k4, k5, k6, k7, k8, v, rates_broadcasts
   use utils
   use constants
   use communication
   use reaclib, only: reaclib_init, reaclib_create_masks, reaclib_preprocessor
   use screening, only: screen_init
   use vital, only: vital_init, vital_rates_derivs
   use netgen, only: iwhatnacre, netgen_init
   use kadonis, only: kadonis_init
   use fuller, only: fuller_init
   use isomers, only: isomers_broadcasts
   use alpha_decays, only: alpha_decays_init, alpha_decays_broadcasts
   use jbj16, only: jbj_init
   use nkk04, only: nkk_init
   use abundances
   use nse_swj, only: nse_init
   use networksetup, only: networksetup_broadcasts
   use bader_deuflhard, only: bader_deuflhard_init
   use evaluate_rates, only: evaluate_all_rates
   use errors, only: MAX_SUBSTEPS
   use reaction_info, only: reaction_info_init
   use physics, only: rnetw2007, rnetw2008
   use reverse, only: reverse_init
   use other_nuc, only: other_nuc_init
   use neutrinos, only: neutrinos_init, neutrino

   implicit none

   real(r8) :: yps(msl,nsp), qi(nre), t9, t9_0, t9_1, t9threshold, ythreshold, &
         rho, rho_0, rho_1, t9max, ye, dzeit, &
         tti1, tti2, tinit, anetw(nsp), znetw(nsp), cpttime, residual
   real(r8), dimension(nsp) :: y0, an, zn
   real(r8), allocatable :: xbuffer(:), xbuffer2(:)

   integer(i4) :: kount, nvar, nvrel, maxnsubt, maxiter, maxtiters, maxnvar1, &
         itime, itraj, pnumrcv, &
         num_sent, isender, stat(MPI_STATUS_SIZE), ierr, &
         num_to_receive, receiver, all_ok, num_fail
   integer(i8) :: ibuffer_size, ibuffer_size2

   character(len=17), parameter :: FMT_meta = "(2I6,ES12.3,4I10)"

   logical :: exit_flag

   type(neutrino) :: nu

   ! ideally replace this common block:
   common / cnetw / an, zn

   ! initialisation
   call tppnp_everything_init()
   call tppnp_output_init()

   if ( master ) then
      call print_column_headers()

      ! get time evolution of tracer particle and its initial abundances; pack and send to slave
      do itraj = 1, num_procs-1
         if ( num_done + num_sent == num_tracers ) exit

         call read_trajectory()
         call set_initial_abundances(yps)
         num_sent = num_sent + 1 ; receiver = num_sent
         call master_pack(yps)
         call master_send(pnum, receiver)
      end do

      ! receive results and send more work
      num_to_receive = num_sent
      do while (num_to_receive > 0)
         call master_receive_and_unpack()
         num_sent = num_sent - 1

         ! check for errors
         select case ( all_ok )
         case(1)
            ! write output
            call write_tppnp_output(pnumrcv, yps, t)
         case(0)
            ! there were errors, do not write
            num_fail = num_fail + 1
            call write_failure(pnumrcv)
         case default
            stop "tppnp: all_ok not understood"
         end select

         ! send more work if there is any
         if (num_done + num_fail + num_to_receive < num_tracers) then
            call read_trajectory()
            call set_initial_abundances(yps)
            num_sent = num_sent + 1 ; receiver = isender
            call master_pack(yps)
            call master_send(pnum, receiver)
            num_to_receive = num_to_receive + 1
         end if
      end do

      ! send exit code to slaves when all work has been sent
      call terminate_all_processes
   end if
   ! master

   if ( slave ) then
      process: do
         call slave_receive_and_unpack(yps)
         if ( exit_flag ) exit ! no more work to do

         call check_mass_conservation(.true., yps, 1, residual)
         if (residual > 1.e-8_r8) call renorm_2d(yps, 1, residual, "")

         ! compute the trajectory
         y0(:) = yps(1,:)
         maxnsubt = 0; maxiter = 0; maxtiters = 0; maxnvar1 = 0
         t%tpeak = ZERO; t%rhopeak = ZERO
         tti1 = wallclocktime()
         do itime = 2, t%ntimesteps
            dzeit = t%time(itime) - t%time(itime - 1)
            t9_0  = t%temp ( itime-1 ) ; t9_1  = t%temp ( itime ) ; t9  = HALF * ( t9_0 + t9_1   ) 
            rho_0 = t%dens ( itime-1 ) ; rho_1 = t%dens ( itime ) ; rho = HALF * ( rho_0 + rho_1 ) 

            !> store peak temperature and density
            t%tpeak = max(t%tpeak, t9_0*1.e9_r8, t9_1*1.e9_r8)
            t%rhopeak = max(t%rhopeak, rho_0, rho_1)

            call calculate_ye( y0, an, zn, ye, considerisotope )

            call evaluate_all_rates( ye, nvar, nvrel, rho, t9, y0, nu )

            call integrate_network( nvar, y0, t9_0, t9_1, rho_0, rho_1, ye, dzeit, nvrel, nu, ierr )

            select case( ierr )
            case( 0 )
               ! all fine
               all_ok = 1
               call check_mass_conservation(.false., y0, residual)
               if (residual > 5.e-7_r8) call renorm_1d(y0, residual, 'integration')
               ! update metadata for time step computation
               maxnsubt  = max(maxnsubt, nsubt)   ; maxiter  = max(maxiter, iter)
               maxtiters = max(maxtiters, titers) ; maxnvar1 = max(maxnvar1, nvar1)
               continue
            case default
               ! an error arose
               all_ok = 0
               exit
            end select

            ! check simulated time, to see if we reached user-defined t_end (defaults to 1e99 s)
            if (t%time(itime) >= t_end) exit
         end do

         yps(1,:) = y0(:)

         cpttime = wallclocktime() - tti1

         ! check status of result
         select case ( all_ok )
         case( 1 )
            ! time integration was good: print tracer particle time integration metadata
            write(*, FMT_meta) pnum, ipid, cpttime, maxnsubt, maxiter, maxtiters, maxnvar1
         case( 0 )
            ! there were issues during time integration, return all_ok = 0 (not okay)
            continue
         end select

         ! send result (i.e. final abundances) to the master
         call slave_pack_and_send_result()

         deallocate(t%time) ; deallocate(t%temp) ; deallocate(t%dens)
      end do process
   end if
   !slave

   if ( master ) write(*,"(a24,es15.6)") "Total program time / s: ", wallclocktime() - tinit
   if ( master ) write(*,*)

   call mpi_barrier(mpi_comm_world, ierr)
   call tppnp_shutdown("Post-processing complete.")

contains

   subroutine tppnp_shutdown(msg)
         character(len=*) :: msg
         if ( master ) then
            !close(ABU_OUTPUT)
            close(TRAJ_INPUT)
         end if
         call shutdown(msg)
   end subroutine tppnp_shutdown

   subroutine slave_pack_and_send_result()
      integer(i4) :: ipos, itag
         itag = pnum
         ipos = 0

         call mpi_pack(all_ok    , 1   , mpi_integer4 , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%mass    , 1   , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%posini  , 3   , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%posfin  , 3   , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%tpeak   , 1   , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%rhopeak , 1   , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)
         call mpi_pack(yps       , nsp , mpi_real8    , xbuffer2 , ibuffer_size2 , ipos , mpi_comm_world , ierr)

         call mpi_send(xbuffer2, ibuffer_size2, mpi_packed, imaster, itag, mpi_comm_world, ierr)
   end subroutine slave_pack_and_send_result

   subroutine slave_receive_and_unpack(yps)
         real(r8) :: yps(nsp)
         integer(i4) :: ipos, ierr

         exit_flag = .false.
         call mpi_recv(xbuffer, ibuffer_size, mpi_packed, imaster, mpi_any_tag, mpi_comm_world, stat, ierr)
         if (stat(mpi_tag) == 0) then
            exit_flag = .true.
            return
         end if

         ipos = 0
         call mpi_unpack(xbuffer , ibuffer_size , ipos , pnum         , 1 , mpi_integer4 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%ntimesteps , 1 , mpi_integer4 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%mass       , 1 , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%posini     , 3 , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%posfin     , 3 , mpi_real8    , mpi_comm_world , ierr)

         if ( allocated( t%time) ) then
            deallocate( t%time ) ; deallocate( t%temp ) ; deallocate( t%dens )
         end if
         allocate( t%time(t%ntimesteps), t%temp(t%ntimesteps), t%dens(t%ntimesteps))

         t%time(:) = ZERO ; t%temp(:) = ZERO ; t%dens(:) = ZERO

         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%time , t%ntimesteps , mpi_real8 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%temp , t%ntimesteps , mpi_real8 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , t%dens , t%ntimesteps , mpi_real8 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer , ibuffer_size , ipos , yps    , nsp          , mpi_real8 , mpi_comm_world , ierr)
   end subroutine slave_receive_and_unpack

   subroutine master_receive_and_unpack()
         integer(i4) :: ipos

         !> receive a message from a slave
         ipos = 0
         call mpi_recv(xbuffer2, ibuffer_size2, mpi_packed, mpi_any_source, mpi_any_tag, mpi_comm_world, stat, ierr)
         num_to_receive = num_to_receive - 1
         isender = stat(mpi_source)

         !> unpack the received message into its constituent parts
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , all_ok    , 1   , mpi_integer4 , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , t%mass    , 1   , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , t%posini  , 3   , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , t%posfin  , 3   , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , t%tpeak   , 1   , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , t%rhopeak , 1   , mpi_real8    , mpi_comm_world , ierr)
         call mpi_unpack(xbuffer2 , ibuffer_size2 , ipos , yps       , nsp , mpi_real8    , mpi_comm_world , ierr)

         pnumrcv = stat(mpi_tag)
   end subroutine master_receive_and_unpack

   subroutine tppnp_output_init()
         implicit none
         integer(i4) :: i, io, io_traj, io_eof, ierr, pos_traj, pos_eof, traj_num
         integer(i8) :: read_pos
         character(len=3) :: eof_string
         logical :: output_file_existed, donefile_existed, failfile_existed

         if ( master ) then

            ! initialise array containing numbers of already-computed particles
            allocate( tracers_done(num_tracers) )
            tracers_done(:) = -1
            num_done = 0

            ! check for a restart
            inquire( file = trim(outfile), exist = output_file_existed )

            ! check whether output file is empty (code ran but no particles computed)
            inquire( file = trim(donefile), exist = donefile_existed )

            if ( output_file_existed .and. donefile_existed) then
               ! this is a restart

               restart = .true.
               print *, 'restart file detected ... '
               print *,

               ! read header of existing output file
               open( file = outfile, unit = ABU_OUTPUT, form = "unformatted", &
                     status = "old", action = "read", access="stream")
               read( ABU_OUTPUT, iostat = io ) nspecies
               if ( io /= 0 ) then
                  write(*,*) "tppnp_output_init: header of existing binary output file wrong, stopping"
                  stop
               endif

               ! check network size in restart file matches current network size
               if ( nspecies /= count(considerisotope) ) then
                  stop "tppnp_output_init: network changed on restart. This is not allowed."
               end if

               ! skip over A and Z arrays
               inquire(ABU_OUTPUT, pos=read_pos)
               read_pos = read_pos + (2 * r8size * nspecies)

               ! loop over number of tracer particles in existing output file reading
               ! the trajectory number, skipping over the result and checking the 3-character
               ! end-of-entry flag is intact
               do i = 1, num_tracers

                  ! read tracer particle number
                  read( ABU_OUTPUT , iostat = io_traj, pos = read_pos ) traj_num
                  inquire( ABU_OUTPUT, pos = read_pos )

                  ! skip over metadata (mass, posini, posfin, tpeak, rhopeak) and
                  ! abundances by seeking to where end of particle marker should be
                  read_pos = read_pos + r8size*(1 + 3 + 3 + 1 + 1 + nspecies)
                  read( ABU_OUTPUT , pos=read_pos, iostat = io_eof  ) eof_string

                  ! check tracer particle was good
                  if ( (io_traj == 0) .and. (io_eof == 0) .and. (eof_string == "EOF") ) then
                     ! particle is good, so count it as already done and store the particle number
                     num_done = num_done + 1
                     tracers_done(i) = traj_num
                     ! save last offset for starting output writing if the next particle is messed up
                     inquire(ABU_OUTPUT,pos=read_pos)
                     last_write_offset = read_pos
                  else
                     ! something is wrong with the tracer particle, stop here and start computing the
                     ! remaining particles
                     exit
                  endif

               end do

               ! close the output file
               close( ABU_OUTPUT )

               write(*,*) "   Particles completed: ", num_done
               write(*,*) "   Resume writing from: ", last_write_offset
            else
               ! this is a fresh run
               restart = .false.

               print *, 'No restart file detected ... '
               nspecies = count(considerisotope)

               ! write network information into output file header
               open( file = outfile, unit = ABU_OUTPUT, form = "unformatted", action = "write",  access="stream")

               write( ABU_OUTPUT ) count(considerisotope) 
               write( ABU_OUTPUT ) pack( an, mask=considerisotope ) 
               write( ABU_OUTPUT ) pack( zn, mask=considerisotope ) 

               ! save last offset for writing
               inquire(ABU_OUTPUT, pos = last_write_offset)
               close( ABU_OUTPUT )
            end if
         end if

         ! delete the failures file if it existed
         inquire( file = trim(failfile), exist = failfile_existed )
         open(file = failfile, unit = FAIL_OUTPUT, status="old", iostat=ierr)
         if (ierr == 0) close(FAIL_OUTPUT, status='delete')

         call broadcast( num_done ) ; call broadcast( num_tracers )

         ! check if everything has already been computed
         if ( num_tracers == num_done ) call tppnp_shutdown( "post-processing complete" )
   end subroutine tppnp_output_init

   subroutine tppnp_init_comm_size()
         num_sent = 0
         call broadcast( nsteps_max )
         ! allocate send buffer size of master->slave communication according to max. steps
         ! mass(1), posini(3), posfin(3), yps(nsp),
         ! (time(nsteps), temp(nsteps), dens(nsteps)), pnum, ntimesteps
         ibuffer_size = (1 + 2*3 + nsp + (nsteps_max * 3)) * 8 + 4 + 4

         ! slave->master communication buffer size
         ! all_ok, mass(1), posini(3), posfin(3), tpeak(1), rhopeak(1), yps(nsp)
         ibuffer_size2 = 4 + (1 + 3 + 3 + 1 + 1 + nsp) * 8
         allocate(xbuffer(ibuffer_size))
         allocate(xbuffer2(ibuffer_size2))
   end subroutine tppnp_init_comm_size

   subroutine write_tppnp_output( pnumrcv, yps, t )
         use trajectories
         ! itraj: particle number
         type(trajectory) :: t
         integer(i4) :: pnumrcv, fh
         real(r8) :: yps(nsp)

         ! open output file
         open( file = outfile, unit = ABU_OUTPUT, form = "unformatted", &
               status = "old", action = "write" , ACCESS="STREAM")   

         ! write the result to the file in the correct position
         write( ABU_OUTPUT , pos = last_write_offset) pnumrcv
         write( ABU_OUTPUT ) t%mass
         write( ABU_OUTPUT ) t%posini(:)
         write( ABU_OUTPUT ) t%posfin(:)
         write( ABU_OUTPUT ) t%tpeak
         write( ABU_OUTPUT ) t%rhopeak
         write( ABU_OUTPUT ) pack( yps, mask=considerisotope )
         write( ABU_OUTPUT ) "EOF"

         ! save last offset for writing
         inquire(ABU_OUTPUT,pos=last_write_offset)

         ! close the output file
         close( ABU_OUTPUT )

         ! count the particle as being done
         num_done = num_done + 1

         ! write num_done to ascii file
         fh = get_unused_file_handle()
         open( file = donefile, unit = fh, action = 'write')
         write( fh, * ) num_done
         close( fh )
   end subroutine write_tppnp_output 

   subroutine write_failure(pnumrcv)
         use traj_data, only: FAIL_OUTPUT, failfile
         integer(i4) :: pnumrcv
   
         open( file = failfile, unit = FAIL_OUTPUT, action = "write", &
               position = "append" )

         write( FAIL_OUTPUT, "(i9)" ) pnumrcv

         close( FAIL_OUTPUT )
   end subroutine write_failure

   subroutine read_input_decks()
         call readframeinput(t9threshold, ythreshold, datdir)
         call readsolverinput()
         call readphysicsinput()
   end subroutine read_input_decks

   subroutine tppnp_everything_init()
         use networksetup, only: read_networksetup, write_networksetup, write_networksetup2, nvnc1, nrnc1
         use array_sizes, only: nsp, i325dim, iCfdim, i282dim, iAtdim
         use nuc_data, only: niso
         tinit = wallclocktime()
         call comm_init()

         if (num_procs == 1) stop "must use more than 1 processor"

#if pIDX_RCLB == 3
         allocate(niso(0:i325dim,0:iCfdim,2))
#else
         allocate(niso(0:i282dim,0:iAtdim,2))
#endif

         call read_input_decks()
         if (master) then
            if (do_neutrinos) stop "neutrinos not fully supported in tppnp yet, sorry"
            if (t_end /= 1.e99_r8) then
               write(*,"(a23,es12.3,a8)") "simulation will end at ", t_end, "seconds"
            end if
         end if

         ! initialise some physics and solver modules
         call reaction_info_init() ; call solver_init()  ; call screen_init()
         call vital_init()         ; call kadonis_init() ; call fuller_init()
         call reaclib_init()       ; call netgen_init()  ; call rates_init()
         call alpha_decays_init()  ; call jbj_init()     ; call nkk_init()
         call bader_deuflhard_init(); call other_nuc_init()
         call neutrinos_init(nu)

         ! set up the network libraries
         istart = 0

         select case(ininet)
         case(3)
            call read_networksetup(an, zn, t9_nw_ini, rho_nw_ini)
            call rnetw2008(ye_nw_ini, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini, nu)
         case(4)
            ! pure reaclib
            stop 'rnetw_reaclib (ininet=4) in progress'
         case default
            call rnetw2007(ye_nw_ini, qi, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini)
         end select
         istart = 1

         if (master) then
            call opening_statement()
            datdir = trim(datdir) // '/'
            write(*,*) "trajectory files are in: ", trim(datdir)
            write(*,*) 'Input files read ...'
            write(*,*)
         end if

         ! build the map from dense to sparse matrix formats
         call makemap()

         ! anetw, znetw are the mass numbers and proton numbers
         ! in yps. an, zn are the
         ! same but for every isotope in networksetup.txt, even the false ones.
         anetw(:) = 0; znetw(:) = 0
         kount = count(considerisotope(:))
         anetw    ( 1 : kount )  = pack ( an             ( 1 : nsp ) , considerisotope ( 1 : nsp )  ) 
         znetw    ( 1 : kount )  = pack ( zn             ( 1 : nsp ) , considerisotope ( 1 : nsp )  ) 

         ! initialise nuclear data
         call nuc_data_init()
         if (detailed_balance) call reverse_init()

         ! ^_^ write networksetup file
         select case(ininet)
         case(3)
            ! write networksetup2
            stop "ppn frame: check souce code arguments to write_networksetup2"
         case(1)
            call write_networksetup(t9, rho, nvnc1, nrcp, nrnc1, an, zn)
            print *, 'stopping in rnetw2007:'
            print *, 'the file networksetup.txt has been created.'
            print *, 'you can make changes to the network by'
            print *, 'modifying this file. When you are happy with'
            print *, 'your changes, run the code again with ininet=3'
            print *, 'and your new networksetup.txt will be used.'
            stop
         case default
            call write_networksetup(t9_nw_ini, rho_nw_ini, nvnc1, nrcp, nrnc1, an, zn)
         end select

         ! create reaclib masks
         call reaclib_create_masks()
         if ( use_cache ) call reaclib_preprocessor

         call nse_init()

         if ( master ) then
            call print_header()
            call traj_data_init()
            num_fail = 0
         end if

         ierr = 0
         if ( master .and. num_procs > num_tracers + 1 ) ierr = 1
         call allreduce_max( ierr )
         if ( ierr /= 0 ) call tppnp_shutdown("too many processors")

         call tppnp_init_comm_size()
   end subroutine tppnp_everything_init


   subroutine terminate_all_processes
         integer :: i
         do i = 1, num_procs - 1
            call terminate_process(i)
         end do
   end subroutine terminate_all_processes


   subroutine master_pack(yps)
         ! pack a message to be sent to a slave
         real(r8) :: yps(nsp)
         integer(i4) :: padsize, ierr, ipos
         integer(i4), allocatable :: packer(:)

         ipos = 0
         call mpi_pack(pnum         , 1            , mpi_integer4 , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%ntimesteps , 1            , mpi_integer4 , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%mass       , 1            , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%posini     , 3            , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%posfin     , 3            , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%time       , t%ntimesteps , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%temp       , t%ntimesteps , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(t%dens       , t%ntimesteps , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)
         call mpi_pack(yps          , nsp          , mpi_real8    , xbuffer , ibuffer_size , ipos , mpi_comm_world , ierr)

         ! ^_^ pad message with zeros
         padsize = int(ibuffer_size - ipos,i4) / 4
         if (padsize > 0) then
            allocate(packer(padsize))
            packer(:) = 0
            call mpi_pack(packer, padsize, mpi_integer4, xbuffer, ibuffer_size, ipos, mpi_comm_world, ierr)
            deallocate(packer)
         end if
   end subroutine master_pack

   subroutine master_send(pnum, receiver)
         integer(i4) :: pnum, receiver
         call mpi_send(xbuffer, ibuffer_size, mpi_packed, receiver, pnum, mpi_comm_world, ierr)
   end subroutine master_send

   subroutine print_header()
         implicit none

         print *, ''
         print *, ''
         print *, "       88P'888'Y88 888 88e  888 88e  Y88b Y88 888 88e  "
         print *, "       P'  888  'Y 888 888D 888 888D  Y88b Y8 888 888D "
         print *, "           888     888 88^  888 88^  b Y88b Y 888 88^  "
         print *, "           888     888      888      8b Y88b  888      "
         print *, "           888     888      888      88b Y88b 888      "
         print *, ''
         print *, ''
         print *, ""
         print *, "       Definitions:"
         print *, "       -----------"
         print *, "       pnum             particle/trajectory number"
         print *, "       ipid             process number"
         print *, "       cpttime          computation time for trajectory"
         print *, "       maxnsubt         max. subtimesteps encoutered"
         print *, "       maxiter          max. iters (per successful [sub]step)"
         print *, "       maxtiters        max. tot. (incl. substeps) iters"
         print *, "       maxnvar1         max. No. species in adaptive network"
         write(*, *)
         write(*, "(A1)", advance = 'no') ' '
         write(*,*)
   end subroutine print_header 


   subroutine print_column_headers
         implicit none
         character*256 :: FMT_head
         integer :: i
         write(*, *)
         FMT_head = "(A6,A6,A12,A10,A10,A10,A10)"
         write(*, FMT_head) 'pnum', 'ipid', 'cpttime', 'maxnsubt', 'maxiter', 'maxtiter', 'maxnvar1'
         do i = 1, 64
            write(*, "(A1)", advance = 'no') '-'
         end do
         write(*, *)

   end subroutine print_column_headers


   subroutine opening_statement()

         write(*,*) "This is TPPNP: trajectory post-processing network."
         write(*,*) " (c) 2016 Marco Pignatari, Falk Herwig,"
         write(*,*) "     Michael Bennett, Raphael Hirschi and the"     
         write(*,*) "     NuGrid collaboration."
         write(*,*)

   end subroutine opening_statement

end program

