program mppnp

   ! *** MPPNP: this is the multi-zone post-processing network/parallel 
   !     frame, which is part of the PPN code suite developed and maintained
   !     by the NuGrid collaboration (http://forum.astro.keele.ac.uk:8080/nugrid)
   !     see user documentation in ../DOC
   ! *** (c) Marco Pignatari, Falk Herwig, Michael Bennett 
   !         and Raphael Hirschi (2007 - 2009), using the USEEPP library
   !         developed by Steven Diehl and Gabriel Rockefeller
   !
   !     See end of this file (mppnp.f90) for variable definitions

   use array_sizes
   use frame_knobs
   use nuc_data
   use physics_knobs
   use solver_knobs
   use solver, only: solver_broadcasts, solver_init, integrate_network
   use solver_diagnostics
   use nse_solver
   use backward_euler
   use jac_rhs
   use rates, only: rates_init, v, rates_broadcasts
   use utils
   use communication
   use screening, only: screen_init
   use reaclib, only: reaclib_init, reaclib_create_masks, reaclib_preprocessor
   use vital, only: vital_init
   use netgen, only: iwhatnacre, netgen_init
   use kadonis, only: kadonis_init
   use fuller, only: fuller_init
   use mixing, only: sig_term_limit, do_mixing, mixing_time
   use fse_wrapper
   use isomers, only: isomers_broadcasts
   use constants
   use alpha_decays, only: alpha_decays_init, alpha_decays_broadcasts
   use jbj16, only: jbj_init
   use nkk04, only: nkk_init
   use nse_swj, only: nse_init
   use bader_deuflhard, only: bader_deuflhard_init
   use evaluate_rates, only: evaluate_all_rates
   use errors, only: MAX_SUBSTEPS, SMALL_DT
   use reaction_info, only: reaction_info_init
   use physics, only: rnetw2007, rnetw2008
   use networksetup, only: read_networksetup, nvnc1, nrnc1, write_networksetup
   use reverse, only: reverse_init
   use other_nuc, only: other_nuc_init
   use neutrinos, only: neutrinos_init, neutrino
   implicit none

   integer, parameter :: maxmod = 1000000, maxfile = 10000, idxfile_fh = 8, elemdim = 120, dimens = 100
   integer :: m, m_old, nvar, nvrel, idum, i, j, icountmod, ipacket, mi, k, idx, ii, iii, jj, mwd, &
         kount, maxzonselt, nprno, nprnr, modstart, modstop, irestart, igrid, idatastepsize, fid1, &
         fid2, ifilemax, nprnlast, ipointofilenr(maxmod), ifirstcycle(maxfile), &
         ilastcycle, isurf, num_of_elements, num_of_elem_print, &
         ierr
   integer, dimension(msl) :: aiter, anvar1, ansubt, atiters
   integer, allocatable :: elem_index(:,:), elem_num(:), elem_index_print(:) 

   real(r8) :: yps(msl,nsp), qi(nre), yps_old(msl,nsp), xmwd, dzeit, dzeitj, ye, &
         t9, t9_0, t9_1, rho, rho_0, rho_1, t9threshold, ythreshold, &
         tti1, tti2, age, xmrmax, xmrmaxi, xmrmin, xfind, dconst, dxm, age_unit, mass_unit, &
         radius_unit, rho_unit, temperature_unit, dcoeff_unit, mini, zini, rotini, overini, ylin, &
         xfac, residual, max_residual, amr_time, rtime, nuctime, surf_coord, xmtest(3)
   real(r8), dimension(nsp) :: y0, an, zn, y_accrete, yps_av
   real(r8), dimension(msl) :: dse, rhose, t9se, rse, xm, dppg, rhoppg, t9ppg, rppg, xmm, &
         xm_old, dq, dq_reverse, antime, artime, anuctime, as
   real(r8), dimension(elemdim) :: elem_surf_decay, xnum_surf_decay, elem_surf_nodecay, xnum_surf_nodecay

   logical :: write_restart, zone_to_compute(msl), file_exists

   character ( len = 30  ) :: filenam(maxfile)
   character ( len = 256 ) :: datdir, filein, fileout(3), prefix, ini_file_name2
   character ( len = 7   ) :: CMODELL
   character ( len = 11  ) :: modname
   character ( len = 80  ) :: codev

   ! ideally remove this common block:
   common / cnetw    / an, zn


   integer(kind=8), parameter :: buffer_size = (nsp + 2)*8, buffer_size2 = nsp*8      
   real(r8) :: xbuffer(buffer_size), xbuffer2(buffer_size2)
   integer :: num_sent, iposition, isender, stat(mpi_status_size), tag, imass_shell, irecv
   logical lsentflg

   ! *** sub-time steping and outputing
   !     isubmax          number of sub-time steps
   !     ksubc            counter variable for sub-time steps
   !     cksubc           name extension character for sub-time step output
   !     icountwsub       counts all timesteps, incl. time-steps and sub-time steps 
   !     writeoutsubcycle if 1, write sub-time steps into all output
   !     icountmodw       similar to icountmod, but counts also sub-time
   !     		       steps if writeoutsubcycle is set to 1
   !     writecycle       true: write output else: skip sub-timestep output 
   ! * in order to activate sub-time steps with output set isubmax to any value >1
   ! * then, also make sure that the statement to create sksubc has the correct format,
   !   e.g. "write(cksubc,'(I2.2)')ksubc" for 2 digits corresponding to isubmax=99
   ! ***********************************************************************************

   integer ksubc, isubmax, icountwsub, writeoutsubcycle, icountmodw
   logical writecycle
   character*4 cksubc
   !HIF (Hydrogen Ingestion Flash )
   integer ing, ingperiod, hif
   real(r8) :: xingest, xHe3Ingest

   ! *** restart from timing parameters
   !     ttstart          time at last restart output (or initialisation)
   !     trestart         duration between restart outputs (input; it determines the
   !                      frequency of the restart output)
   !     tstop            time at which the program should stop (input)
   ! ***********************************************************************************
   integer :: ioutc, ioutformat, isomnetw(nsp)
   real(r8) :: tinit, tfin, ttstart, tcurrent, tstop, ttdiff, trestart, anetw(nsp), znetw(nsp)
   real(r8), allocatable :: yps_restart(:,:), yps_out(:,:)
   real(r8), allocatable, dimension(:) :: yps_out_decay, yps_out_nodecay, &
         elem_out_decay, xnum_out_decay, elem_out_nodecay, xnum_out_nodecay
   ! if terminate_flag == .true., the program terminates after the timestep
   logical :: terminate_flag
   ! v_terr are the terrestrial rates used for subroutine decay, to
   ! calculate in the envelope decayed abundances (yps_surf_decay).
   real(r8) :: yps_surf_decay(nsp), yps_surf_nodecay(nsp), v_terr(nre), t9_for_decay

   ! neutrino object:
   type(neutrino) :: nu

   call comm_init

   if ( master ) then
      tinit = wallclocktime() ; ttstart = wallclocktime()
      call opening_statement()
   end if

   if ( num_procs ==  1 ) stop "At least two processors are required. Please increase from 1"

#if pIDX_RCLB == 3
   allocate(niso(0:i325dim,0:iCfdim,2))
#else
   allocate(niso(0:i282dim,0:iAtdim,2))
#endif

   call readframeinput(t9threshold, ythreshold, modstart, modstop, igrid, &
         dxm, xmrmin, xmrmax, datdir, prefix, nprno, nprnr, &
         ioutc, trestart, tstop, hif, sig_term_limit, isubmax, writeoutsubcycle)
   call readsolverinput()
   call readphysicsinput()

   ! initialisations without dependencies on input files
   call reaction_info_init() ; call solver_init()  ; call bader_deuflhard_init()
   call vital_init()         ; call kadonis_init() ; call fuller_init()
   call reaclib_init()       ; call netgen_init()  ; call rates_init()
   call alpha_decays_init()  ; call jbj_init()     ; call nkk_init()
   call other_nuc_init()     ; call neutrinos_init(nu)

   dq(:) = ZERO ; isurf = 0 ; v_terr(:) = ZERO ; terminate_flag = .false.


   if (hif == 1) call readhifinput(ing, ingperiod, xingest)

   if ( master ) then

      if (do_neutrinos) stop "neutrinos not fully integrated into mppnp yet, sorry"
      open(summary_output, file = "summaryinfo.dat")

      ! setting up file names/lists... to be moved
      ! TODO: all of the following code that concerns input file names should probably move to a
      ! module or separate routine that deals with that
      datdir = trim(datdir) // '/'
      if (iolevel >= 3) print *,"USEEPP input files are in: ", trim(datdir) // 'end'

      ! USEEPP cyclenb/file pointer array set-up
      filein = trim(datdir) // trim(prefix) // ".idx"
      call check_file_exists(filein)
      open(idxfile_fh, file = filein, action = 'read' )
      i = 0
      do while (.true.)
         i = i + 1
         read(idxfile_fh,*,end=9833) filenam(i)
      end do
      9833 close(idxfile_fh)
      ifilemax = i - 1

      if (iolevel >= 2) print *, "There are ", ifilemax, " h5 packet files to be processed."

      ! read model numbers from filenames in *.idx 
      do i = 1, ifilemax
         read(filenam(i)(int(len_trim(prefix))+2:int(len_trim(prefix))+8),'(I7)') ifirstcycle(i)
      enddo

      ! to be destroyed by fire and done properly:
      ! because FH messed up writing the right number into icyclenb of his
      ! stellar evolution output (he wrote nprnlast into icyclenb,
      ! how stupid is that) we basically can not use the last cycle). for
      ! now we don't have to. one workaround would be to add a fake line
      ! to the idx file. but for now just don't use the data in the last
      ! file and take:
      nprnlast = ifirstcycle(ifilemax)-1
      ilastcycle = nprnlast
      if (modstop > 0) ilastcycle = min(modstop,nprnlast)
      if (iolevel >= 2) print *, "First and last cycle of available input models: ", ifirstcycle(1), nprnlast

      ! now setup pointer array to find the h5 file in which a particular cycle is located       
      do i = 1, ifilemax-1
         ii  = ifirstcycle(i)
         iii = ifirstcycle(i+1) - 1
         if (iii > maxmod) stop "MPPNP: hardcoded parameter maxmod too small"
         ipointofilenr(ii:iii) = ifirstcycle(i)
      end do
      ii  = ifirstcycle(ifilemax) 
      iii = nprnlast
      ipointofilenr(ii:iii) = ifirstcycle(ifilemax)

      if (iolevel >= 4) then
         do i = modstart,iii
            print *, "i, ipointofilenr(i) ", i, ipointofilenr(i)
         end do
      endif

      ! end of file name reading and list making...

      write(*,*) 'Input files read ...'
      write(*,*) 'setting up network ...'
   end if ! master

   istart = 0
   select case(ininet)
   case(3)
      call read_networksetup(an, zn, t9_nw_ini, rho_nw_ini)
      call rnetw2008( ye_nw_ini, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini, nu )
   case default
      call rnetw2007( ye_nw_ini, qi, an, zn, nvar, nvrel, rho_nw_ini, t9_nw_ini, yps_nw_ini )
   end select
   istart = 1

   ! make mapping from dense jacobian to CSR (row-compressed) jacobian
   call makemap()

   ! deal with restarting: also needs to be moved to its own routine out of the way
   if (master) then
      ! *** is this a restart?
      open( restart_log, file = 'last_restart.out' )
      !CR :for sub-time step output, read number of time (incl. sub-time) steps
      if ((isubmax  >  1) .and. (writeoutsubcycle  ==  1)) then
         read(restart_log, *, iostat=ierr) icountmodw, irestart, ipacket
         if  (ipacket  ==  0) then
            print *,"STOP: With sub-time steps output MPPNP requires 3 entries in last_restart.out:"
            print *, "time-step (incl. sub-time step) cycle (1)," // &
                  "SE time-step cycle (2) and first cycle number in restart file (3)." // &
                  "A new sub-time step calculation might then require (1) to be equal (2)."  
            stop
         endif
         if (iolevel >= 3) then
            print *, "Read last_restart.out:", icountmodw, irestart, ipacket
         end if
      else
         read(restart_log,*) irestart, ipacket
      end if
      close(restart_log)
      if (irestart /= 0) then
         iabuini = 0
         modstart = irestart + 1
         write(CMODELL,'(I7.7)') ipacket
         fileout(1) = "H5_restart/" // trim(prefix) // '.' // CMODELL // ".restart.h5"
         fileout(2) = "H5_surf/"    // trim(prefix) // '.' // CMODELL // ".surf.h5"
         fileout(3) = "H5_out/"     // trim(prefix) // '.' // CMODELL // ".out.h5"
      end if

      ! *** read in 1D stellar evolution output for post-processing ************
      write(filein,'(i7.7)') ipointofilenr(modstart)
      filein = trim(datdir) // trim(prefix) // "." // trim(filein) // ".se.h5"
      call check_file_exists( filein )

      call FSE_OPEN( filein, FID1 )

      ! *** reading in global parameters:
      !     Units for variables are with respect to CGS
      call FSE_READ_SATTR( FID1 , -1 , "codev"            , codev            ) 
      call FSE_READ_SATTR( FID1 , -1 , "modname"          , modname          ) 
      call FSE_READ_DATTR( FID1 , -1 , "mini"             , mini             ) 
      call FSE_READ_DATTR( FID1 , -1 , "zini"             , zini             ) 
      call FSE_READ_DATTR( FID1 , -1 , "rotini"           , rotini           ) 
      call FSE_READ_DATTR( FID1 , -1 , "overini"          , overini          ) 

      call FSE_READ_DATTR( FID1 , -1 , "age_unit"         , age_unit         ) 
      call FSE_READ_DATTR( FID1 , -1 , "mass_unit"        , mass_unit        ) 
      call FSE_READ_DATTR( FID1 , -1 , "radius_unit"      , radius_unit      ) 
      call FSE_READ_DATTR( FID1 , -1 , "rho_unit"         , rho_unit         ) 
      call FSE_READ_DATTR( FID1 , -1 , "temperature_unit" , temperature_unit ) 
      call FSE_READ_DATTR( FID1 , -1 , "dcoeff_unit"      , dcoeff_unit      ) 

      if (iolevel  >=  1) write(*,*) "There are ",ilastcycle-modstart, &
            " models to be processed. Now building pp-grid ..."

      ! xmrmaxi = stellar surface of the first model
      call FSE_READ_DATTR ( FID1 , modstart , "total_mass" , xmrmaxi      ) 
      call FSE_READ_IATTR ( FID1 , modstart , "shellnb"    , mi           ) 

      call FSE_READ_D     ( FID1 , modstart , mi           , "mass"       , xmm ) ! Msun
      call FSE_READ_D     ( FID1 , modstart , mi           , "delta_mass" , dq  ) ! Msun

      if (iabuini /= 0) then ! restart
         ! select and initialise a grid for the post-processing
         select case(igrid)
         case(1)
            xmrmax = min(xmrmaxi,xmrmax)
            call customgrid(xmrmin,xmrmax,dxm,xm,m)
         case(2)
            ! using stellar evolution grid       
            ! the nova calculations are only done for a surface layer, from xmrmin.
            ! the first case is accreting of solar composition according to the
            ! file specified according to iabuini
            if ( xmrmin < ZERO .and. iabuini == 10 ) xmrmin = 0.999_r8 * xmm(1)
            if ( xmrmin < ZERO .and. iabuini == 20 ) then
               xmwd   = -xmrmin
               xmrmin = 0.999_r8 * xmwd
            endif
            if (iabuini == 20) then  ! the following lines are about! fining mwd?!
               i=1
               do while ( xmm(i) > xmwd.and.i < mi )
                  i=i+1
               enddo
               mwd=i
            end if
            i=1          ! next below find index of xmrmin in xmm
            do while ( xmm(i) > xmrmin.and.i < mi )
               i=i+1
            enddo
            mi=i
            do i=1,mi
               dq_reverse(i)=dq(mi+1-i)
               xm(i)=xmm(mi+1-i)
            enddo
            m=mi
            ! *** dpa ***  we reverse indexing, so need to update mwd
            if (iabuini == 20) mwd = m+1-mwd
         case(3)
            stop 'amr depreciated 2016-09-25'
         end select
      end if
      ! for subtimestep output, restart from cycle number icoundmow
      if ( (isubmax > 1) .and. (writeoutsubcycle == 1)) irestart = icountmodw

      ! set the initial abundances for the whole star
      call iniabund(m, xm ,y_accrete, mwd ,yps, irestart, fileout(1), nvar)

   end if ! master

   ! anetw, znetw and isomnetw are the A, Z and isomeric states of the isotopes that are true in
   ! networksetup.txt
   anetw(:) = ZERO ; znetw(:) = ZERO ; isomnetw(:) = 0

   kount = count(considerisotope(1:nsp))
   anetw    ( 1:kount )  = pack ( an             ( 1:nsp )  , considerisotope ( 1:nsp )  ) 
   znetw    ( 1:kount )  = pack ( zn             ( 1:nsp )  , considerisotope ( 1:nsp )  ) 
   isomnetw ( 1:kount )  = pack ( isomeric_state ( 1:nsp )  , considerisotope ( 1:nsp )  ) 

   num_of_elements = nint(maxval(znetw)) + 1

   if (master) then
      ! for surf_elem routine
      allocate ( yps_out_nodecay  ( nvar            )        ) 
      allocate ( yps_out_decay    ( nvar            )        ) 
      allocate ( elem_num         ( num_of_elements )        ) 
      allocate ( elem_index_print ( num_of_elements )        ) 
      allocate ( elem_out_nodecay ( num_of_elements )        ) 
      allocate ( xnum_out_nodecay ( num_of_elements )        ) 
      allocate ( elem_out_decay   ( num_of_elements )        ) 
      allocate ( xnum_out_decay   ( num_of_elements )        ) 
      allocate ( elem_index       ( num_of_elements , dimens )  ) 
   end if 

   tti1 = wallclocktime()

   call broadcast(modstart)   ;  call broadcast(ilastcycle)

   ! initialise modules that do require that the network is initialized
   call nuc_data_init()
   call screen_init()
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

   ! this is the cycle loop

   if (master .and. iolevel >= 3) write(*,*) "Cycle loop starts with modstart, ilastcycle= ",modstart, ilastcycle

   ! set up icountwsub counter
   icountwsub = 0
   if ( (isubmax  >  1) .and. (writeoutsubcycle  ==  1) ) then
      if (irestart == 0) then
         icountwsub = 0
      else
         !in case of a restart
         icountwsub = icountmodw - modstart +1
      end if
   end if

   cycleloop: do icountmod = modstart, ilastcycle

      if ((isubmax  >  1) .and. (writeoutsubcycle  ==  1)) then
         icountwsub = icountwsub -1
      end if

      subcycle: do ksubc = 0, (isubmax-1)

         !CR: for sub-time step parameter
         writecycle = .true.
         icountmodw = icountmod
         if  ( isubmax  >  1 ) then
            if ( writeoutsubcycle  ==  1 ) then
               icountwsub = icountwsub + 1
               icountmodw = icountmod + icountwsub
            else
               writecycle = .false.
               if ( ksubc  ==  (isubmax-1) ) writecycle = .true.
            end if
         end if

         max_residual = ZERO

         ! Now set master part.  Slave part skips to the broadcasts below
         if (master) then
            tti1 = wallclocktime()
            !CR: for sub-time step output
            if (isubmax > 1) then
               write(summary_output,'(A15,I6,A15,I6,A26,ES12.2)') &
                     "# Next cycle = ", icountmod,", subtimestep = ",ksubc, &
                     ", run time to this point: ",tti1
            else
               write(summary_output,'(A15,I6,A26,ES12.2)') "# Next cycle = ", &
                     icountmod,", run time to this point: ",tti1
            end if

            write(summary_output,'(A59,A25)') "# M | iter | titer | nsubt | < NRNW t/s > | nucnet99 t/s | ", &
                  "last nvar1 | physics t/s "

            !     open file for h5 output: 
            idum = 0
            i = icountmodw - ipacket + ioutc
            if (modulo(i,ioutc) == 0 .or. (irestart == 0 .and. icountmodw == modstart)) then
               ipacket = icountmodw
               idum = 1
            else if (irestart > 0 .and. iabuini == 0) then
               idum = 1
               irestart = 0
            end if

            if (idum == 1 .and. writecycle) then

               write(CMODELL,'(I7.7)') ipacket

               ! The following output files are used in MPPNP.  They are written every
               ! ioutc cycles.  The header information is the same in each case, so we
               ! can just use a loop
               fileout(1) = "H5_restart/" // trim(prefix) // '.' // CMODELL // ".restart.h5"
               fileout(2) = "H5_surf/" // trim(prefix) // '.' // CMODELL // ".surf.h5"
               fileout(3) = "H5_out/" // trim(prefix) // '.' // CMODELL // ".out.h5"
               do j = 1, 3
                  call FSE_OPEN(fileout(j),FID2)
                  ! *** writing global parameters:
                  call FSE_WRITE_IATTR(FID2, -1, "numcodev",ioutc    )
                  call FSE_WRITE_SATTR(FID2, -1, "codev",   codev    )
                  call FSE_WRITE_SATTR(FID2, -1, "modname", modname  )
                  call FSE_WRITE_DATTR(FID2, -1, "mini",    mini     )
                  call FSE_WRITE_DATTR(FID2, -1, "zini",    zini     )
                  call FSE_WRITE_DATTR(FID2, -1, "rotini",  rotini   )
                  call FSE_WRITE_DATTR(FID2, -1, "overini", overini  )
                  call FSE_WRITE_IATTR(FID2, -1, "zisnb",   nvar     )

                  ! *** write full output if it is a restart file
                  call FSE_WRITE_DARRAYATTR(FID2, -1, "Z",              znetw,    nvar)
                  call FSE_WRITE_DARRAYATTR(FID2, -1, "A",              anetw,    nvar)
                  call FSE_WRITE_IARRAYATTR(FID2, -1, "isomeric_state", isomnetw, nvar)
                  ! *** Units for variables (with respect to CGS)
                  call FSE_WRITE_DATTR(FID2, -1, "age_unit",         age_unit    )
                  call FSE_WRITE_DATTR(FID2, -1, "mass_unit",        mass_unit   )
                  call FSE_WRITE_DATTR(FID2, -1, "radius_unit",      radius_unit )
                  call FSE_WRITE_DATTR(FID2, -1, "rho_unit",         rho_unit    )
                  call FSE_WRITE_DATTR(FID2, -1, "temperature_unit", 1.d9        )
                  call FSE_WRITE_DATTR(FID2, -1, "dcoeff_unit",      dcoeff_unit )     
                  call FSE_CLOSE(FID2)
               end do

            end if

            if (iolevel >= 4) print *, " about to read stellar evolution input"

            ! HIF: add ingestion composition
            ! *** this is a little modification for H-ingestion flash hif
            if (hif  ==  1) then
               ! *** initial restart model contains already first ingestion, so +1
               if (modulo((icountmod - modstart), ingperiod) == 0) then
                  if (iolevel >= 2)  print *, "inserting H for HIF test at model ksubc= ", ksubc
                  do j = 1, m
                     if (xm(j) >= xmrmax - 4.e-4_r8) then
                        xHe3Ingest = xingest * 2.92e-5_r8 / 0.7_r8
                        yps(j,ispe("PROT ")) = yps(j,ispe("PROT ")) + xingest
                        yps(j,ispe("HE  3")) = yps(j,ispe("HE  3")) + xHe3Ingest
                        yps(j,ispe("C  12")) = yps(j,ispe("C  12")) - xingest - xHe3Ingest
                     end if
                  end do
               end if
            end if
            ! *** end hif modification

            ! *** open stellar evolution input file and read stellar structure
            write(filein,'(i7.7)') ipointofilenr(icountmod)
            filein = trim(datdir) // trim(prefix) // "." // trim(filein) // ".se.h5"

            if (iolevel >= 3) then
               print *, 'starting model, from file number, file name:'
               print *, icountmod, ipointofilenr(icountmod), trim(filein)
            endif

            call check_file_exists(filein)
            call FSE_OPEN(filein, FID1)
            call FSE_READ_DATTR ( FID1 , icountmod , "age"     , age           ) 
            call FSE_READ_DATTR ( FID1 , icountmod , "deltat"  , dzeit         ) 
            call FSE_READ_IATTR ( FID1 , icountmod , "shellnb" , mi            ) 
            call FSE_READ_D     ( FID1 , icountmod , mi        , "dcoeff"      , dse   ) ! [d] = cm^2/s
            call FSE_READ_D     ( FID1 , icountmod , mi        , "radius"      , rse   ) ! cm
            call FSE_READ_D     ( FID1 , icountmod , mi        , "mass"        , xmm   ) ! Msun
            call FSE_READ_D     ( FID1 , icountmod , mi        , "delta_mass"  , dq    ) ! Msun
            call FSE_READ_D     ( FID1 , icountmod , mi        , "rho"         , rhose ) ! [cgs]
            call FSE_READ_D     ( FID1 , icountmod , mi        , "temperature" , t9se  ) ! [K]
            call FSE_CLOSE(FID1)

            ! units
            rse   ( 1:mi ) = rse   ( 1:mi )  * radius_unit
            rse   ( 1:mi ) = rse   ( 1:mi )  / 6.96e10_r8 ! convert to Rsun units
            rhose ( 1:mi ) = rhose ( 1:mi )  * rho_unit
            t9se  ( 1:mi ) = t9se  ( 1:mi )  * temperature_unit / 1.e9_r8
            dse   ( 1:mi ) = dse   ( 1:mi )  * dcoeff_unit

            ! adjustments
            if (code_source == 'MES') dzeit = dzeit * age_unit ! (dzeit is in seconds for GENEC)
            dzeitj = dzeit / yrs2sec
            if ( code_source  ==  'GNV' .and. maxval(dse(1:mi))  >  1.e14_r8 ) then
               ! ^_^ what?!...
               if (modulo(icountmod,1000) == 1) print *, 'reducing dcoeff GNV by 10000'
               dse(1:mi) = dse(1:mi) / 1.e4_r8
            endif
            ! end adjustments

            if ( iolevel  >=  2 .or. modulo(icountmod,1000) ==  1 ) write(*,*) "cycle, time step /yrs: ", icountmod, dzeitj

            if ( iolevel >= 3) write(*,fmtprofile) "reading data time/s = '", wallclocktime() - tti1

            select case(igrid)
            case(2)
               tti2 = wallclocktime()
               ! update state with the new state just read from the SE file
               xm_old(:) = xm(:)
               m_old = m

               ! count number of mass shells with mass > xmrmin are there in the se profile
               i = 1
               do while ( xmm(i)  >  xmrmin .and. i  <  mi )
                  i = i + 1
               enddo

               mi = i
               do i = 1, mi
                  dq_reverse  (i)  = dq     (mi+1-i)
                  xm          (i)  = xmm    (mi+1-i)
                  rhoppg      (i)  = rhose  (mi+1-i)
                  dppg        (i)  = dse    (mi+1-i)
                  rppg        (i)  = rse    (mi+1-i)
                  t9ppg       (i)  = t9se   (mi+1-i)
               enddo
               m = mi

               ! interpolate abundances yps from the old grid in the new grid just read from the SE file
               yps_old (:,:) = yps  (:,:)
               yps     (:,:) = ZERO

               do j = 1, m
                  xfind = xm(j)
                  if (xfind <= xm_old(m_old)) then
                     idx = minloc(xm_old(2:m_old), mask = xm_old(2:m_old) >= xfind, dim = 1) + 1
                     xfac = (xfind - xm_old(idx-1)) / (xm_old(idx) - xm_old(idx-1))
                     yps(j,:) = yps_old(idx-1,:) + xfac * (yps_old(idx,:) - yps_old(idx-1,:))
                     ! clip negative abundances introduced by extrapolation in central cell(s)
                     yps(j,:) = max(yps(j,:), 1.e-99_r8)
                  else
                     ! mass has increased; must be due to accretion, so assign it the accretion composition
                     yps(j,:) = y_accrete(:)
                  end if
               end do
            case(1)
               ! *** static grid, interpolate stellar structure onto post-processing grid
               rhoppg (:)  = ZERO
               t9ppg  (:)  = ZERO
               dppg   (:)  = ZERO
               rppg   (:)  = ZERO
               do j = 1, m
                  xfind = xm(j)
                  idx = minloc(xmm(1:mi), mask = xmm(1:mi) > xfind, dim = 1)
                  rhoppg ( j )  = ylin ( xfind, xmm ( idx ) , xmm ( idx + 1 ) , rhose ( idx ) , rhose ( idx + 1 )  ) 
                  t9ppg  ( j )  = ylin ( xfind, xmm ( idx ) , xmm ( idx + 1 ) , t9se  ( idx ) ,  t9se ( idx + 1 )  ) 
                  dppg   ( j )  = ylin ( xfind, xmm ( idx ) , xmm ( idx + 1 ) , dse   ( idx ) ,   dse ( idx + 1 )  ) 
                  rppg   ( j )  = ylin ( xfind, xmm ( idx ) , xmm ( idx + 1 ) , rse   ( idx ) ,   rse ( idx + 1 )  ) 
               end do
            end select

            if (m > MSL) then
               print *, "m (number of zones in stellar model) is larger than msl (max number of mppnp zones)."
               print *, "please increase msl in ppn_frame.input and recompile with `make distclean`, `make`."
               print *, "m = ", m, "msl = ", msl
               stop
            end if

            ! *** decide which zones to actually post-process based on temperature
            call zoneselect(t9threshold, t9ppg, m, zone_to_compute, maxzonselt)

            !CR: for sub-time steps adapt timesteps and ages
            if (isubmax  >  1) then
               age   =  age * age_unit - dzeit * (ONE - (ksubc + 1) / dble(isubmax))
               dzeit =  dzeit / dble(isubmax)
            end if
         end if

         call broadcast(dzeit)     ! time step
         call broadcast(icountmod) ! cycle number

         ! ^_^ init diagnostic arrays
         aiter(:) = 0; atiters(:) = 0; ansubt(:) = 0; anvar1(:) = 0
         antime(:) = ZERO; artime(:) = ZERO; anuctime(:) = ZERO

         if (master) then
            tti1 = wallclocktime()
            ! Master cycle loop: pack state info; send to slaves, receive results, send more work until complete
            print *, 'computing cycle', icountmod

            if ( .not. any( zone_to_compute(1:m) ) ) then
                  call terminate_all_slaves
            else
               ! Send each slave process some initial work
               j = 0
               num_sent = 0
               do while (num_sent < num_procs - 1 .and. j <= m)
                  j = j + 1
                  if ( zone_to_compute(j) ) then
                     num_sent = num_sent + 1
                     rho = rhoppg(j) ; t9 = t9ppg(j) ; y0(:) = yps(j,:)
                     ! pack zone info into the mesage buffer and send to slave <num_sent> with tag <j> (mesh point)
                     call master_pack_and_send(rho, t9, y0, xbuffer, buffer_size, num_sent, j)
                  end if
               end do

               ! kill slaves there's no work for
               if (num_procs - 1  >  num_sent) then
                  do j = num_sent + 1, num_procs - 1
                     call terminate_slave(j)
                  end do
               end if

               ! start receiving results and issuing more work
               irecv = num_sent

               RCV: do while (irecv > 0)
                  call master_recv_and_unpack(xbuffer2, buffer_size2, stat, y0)

                  irecv = irecv - 1 ; isender = stat(mpi_source) ; imass_shell = stat(mpi_tag)

                  ! abundances (yps) are updated with y0, the abundances calculated by the solver for this mesh point
                  yps(imass_shell,:) = y0(:)
                  lsentflg = .true.

                  FNS: do while (.true.) ! find next zone_to_compute
                     if (j >= m) exit FNS
                     j = j + 1

                     if (zone_to_compute(j)) then
                        irecv = irecv    + 1 ; num_sent = num_sent + 1
                        rho   = rhoppg(j)    ; t9       = t9ppg(j)     ; y0 = yps(j,:)
                        call master_pack_and_send(rho, t9, y0, xbuffer, buffer_size, isender, j)
                        lsentflg = .false.
                        exit FNS
                     endif
                  end do FNS

                  ! terminate this slave (for this timestep) if there is no more work to assign him/her
                  if (lsentflg) call terminate_slave(isender)
               end do RCV
            end if
            ! end if (.not. any(zone_to_compute))
         else
            !	This is the slave code: receive, unpack, calculate, pack, send
            do while (.true.)
               call slave_recv_and_unpack(rho, t9, y0, xbuffer, buffer_size, tag)
               if (tag == 0) exit ! no more work to do this time step

               tti1 = wallclocktime()
               call calculate_ye(y0, an, zn, ye, considerisotope)
               call evaluate_all_rates( ye, nvar, nvrel, rho, t9, y0, nu )
               rtime = wallclocktime() - tti1

               t9_0 = t9 ; t9_1 = t9 ; rho_0 = rho ; rho_1 = rho
               tti1 = wallclocktime()
               call integrate_network( nvar, y0, t9_0, t9_1, rho_0, rho_1, ye, dzeit, nvrel, nu, ierr )

               select case( ierr )
               case( 0 ) ! all fine
                  continue
               case( MAX_SUBSTEPS )
                  stop "mppnp: max. substeps reached"
               case( SMALL_DT )
                  stop "mppnp: time step underflow"
               case default
                  stop "mppnp: an exception was raised during time integration"
               end select

               call check_mass_conservation(.false., y0, residual)
               if (residual > 5.e-7_r8) call renorm_1d(y0, residual, 'time integration')
               nuctime = wallclocktime() - tti1

               select case( integration_method )
               case(EULER)
                  ntime = nuctime / titers
               case(DEUFL)
                  ntime = nuctime / max(iter, 1)
               end select


               aiter(tag)  = iter  ; atiters(tag) = titers ; ansubt(tag)   = nsubt   ; anvar1(tag) = nvar1
               artime(tag) = rtime ; antime(tag)  = ntime  ; anuctime(tag) = nuctime

               call slave_pack_and_send( y0, xbuffer2, buffer_size2, tag )
            end do
         end if

         ! Diagnostic information is stored
         call reduce_max ( aiter    )  ;  call reduce_max ( atiters )  ;  call reduce_max ( anvar1 ) 
         call reduce_max ( ansubt   )  ;  call reduce_max ( artime  )  ;  call reduce_max ( antime ) 
         call reduce_max ( anuctime ) 

         call mpi_barrier(mpi_comm_world, ierr)

         if (master) then
            ! mixing
            call do_mixing(xm, dq_reverse, rhoppg, rppg, dppg, yps, m, dzeit)

            ! check for deviations from sum(yps) = 1 after mixing step
            call check_mass_conservation(.false., yps, m, residual)

            if (residual >= SUM_TOL) then
               write(*,fmtprofile) 'WARNING: post-diffusion mass conservation residual: ', residual
               call renorm_2d(yps, m, residual, "diffusion")
            end if

            !CR: with sub-time steps but no sub-time step output jump over printing
            if (.not. writecycle) goto 1337

            ! clip
            yps(:,:) = max(yps(:,:), 1.e-99_r8)

            ! decayed/undecayed surface abundances: a network call with large timestep and low
            ! temperature/density is used for the decay
            yps_surf_decay(:)    = ZERO
            yps_surf_nodecay(:)  = yps(m,:)
            t9_for_decay         = t9ppg(m)

            tti1 = wallclocktime()
            call decay_I( t9_for_decay, nvar, age, nvrel, v_terr, yps_surf_nodecay, yps_surf_decay)
            if ( iolevel >= 3 ) write(*,fmtprofile) "decay time/s = ",wallclocktime() - tti1

            ! compute surface decayed/undecayed elemental abundances
            yps_out_nodecay(:) = pack(yps(m,:), considerisotope)
            yps_out_decay(:)   = pack(yps_surf_decay(:), considerisotope)

            tti1 = wallclocktime()
            call surf_elem(isurf, znetw, anetw, yps_out_nodecay, &
                  nvar,num_of_elements,elem_surf_nodecay, &
                  xnum_surf_nodecay, elem_index,elem_num, dimens)
            call surf_elem(isurf, znetw, anetw, yps_out_decay, &
                  nvar,num_of_elements,elem_surf_decay, xnum_surf_decay, &
                  elem_index,elem_num, dimens)
            if ( iolevel >= 3 ) write(*,fmtprofile) "surf elem time/s = ",wallclocktime() - tti1


            ! write output

            if (iolevel >= 4) print*, "starting output .."

            ! TODO: move to separate output routine
            if (modulo(icountmodw,nprno)  ==  0) then
               ! write full output every nprno cycles
               ! *** write all considered isotopes.
               allocate (yps_out(m,nvar))
               do i = 1, m
                  yps_out(i,:) = pack(yps(i,:), considerisotope)
               end do
               call FSE_OPEN(fileout(3),FID2)
               call FSE_WRITE(FID2,icountmodw, m, 6, &
                     xm, "mass", SE_DOUBLE, &
                     rppg, "radius", SE_DOUBLE, &
                     RHOppg, "rho", SE_DOUBLE, &
                     t9ppg, "temperature", SE_DOUBLE, &
                     dppg, "dcoeff", SE_DOUBLE, &
                     yps_out, "iso_massf", SE_DOUBLE_2D, m, nvar)
               call FSE_WRITE_IATTR(FID2,icountmodw, "shellnb", m)
               call FSE_WRITE_DATTR(FID2,icountmodw,"age", age)
               call FSE_WRITE_DATTR(FID2,icountmodw, "deltat", dzeit)
               call FSE_WRITE_DATTR(FID2,icountmodw, "total_mass", xmrmax)
               call FSE_WRITE_IATTR(FID2,icountmodw, "cycle", icountmodw)
               call FSE_CLOSE(FID2)

               deallocate (yps_out)
            end if
            
            ! TODO: move to separate output routine
            ! surface yields (written every timestep )
            elem_out_nodecay(:)  = elem_surf_nodecay(:num_of_elements)
            xnum_out_nodecay(:)  = xnum_surf_nodecay(:num_of_elements)
            elem_out_decay(:)    = elem_surf_decay(:num_of_elements)
            xnum_out_decay(:)    = xnum_surf_decay(:num_of_elements)

            call FSE_OPEN(fileout(2),FID2)
            call FSE_WRITE(FID2,icountmodw, 1, 11, xm(m), "mass", SE_DOUBLE, &
                  rppg(m), "radius", SE_DOUBLE, RHOppg(m), "rho", SE_DOUBLE, &
                  t9ppg(m), "temperature", SE_DOUBLE, dppg(m), "dcoeff", SE_DOUBLE, &
                  yps_out_nodecay, "iso_massf", SE_DOUBLE_2D, 1, nvar, &
                  elem_out_nodecay,"elem_massf",SE_DOUBLE_2D, 1, num_of_elements, &
                  xnum_out_nodecay, "elem_numf",SE_DOUBLE_2D, 1, num_of_elements, &
                  yps_out_decay, "iso_massf_decay", SE_DOUBLE_2D,1,nvar, &
                  elem_out_decay, "elem_massf_decay", SE_DOUBLE_2D, 1, num_of_elements, &
                  xnum_out_decay, "elem_numf_decay", SE_DOUBLE_2D,1,num_of_elements)
            call FSE_WRITE_IATTR(FID2, icountmodw, "model_number",   icountmodw  )
            call FSE_WRITE_IATTR(FID2, icountmodw, "shellnb",        1           )
            call FSE_WRITE_DATTR(FID2, icountmodw, "age",            age         )
            call FSE_WRITE_DATTR(FID2, icountmodw, "deltat",         dzeit       )
            call FSE_CLOSE(FID2)

            if (iolevel >= 4) print *, "starting output USEEPP .."
            ! USEEPP output: time-dependent restart
            write_restart = .false.
            call cpu_time(tcurrent)
            ttdiff = (tcurrent - ttstart) / 3600._r8
            if (ttdiff > trestart) then
               ttstart = tcurrent
               write_restart = .true.
               if (iolevel >= 1) print *, ttdiff, " hrs passed -> print restart cycle", icountmod
            else if (modulo(icountmodw,nprnr) == 0) then
               write_restart = .true.
            end if

            ! TODO: move to separate output routine
            ! writing restart structure every nprnr models or trestart hours of computed time
            if (write_restart) then
               ! prepare array yps_restart:
               allocate (yps_restart(m,nvar))
               do i = 1, m
                  yps_restart(i,:) = pack(yps(i,:), considerisotope)
               end do

               call FSE_OPEN(fileout(1),FID2)
               call FSE_WRITE(FID2,icountmodw, M, 6, xm, "mass", SE_DOUBLE, &
                     rppg, "radius", SE_DOUBLE, RHOppg, "rho", SE_DOUBLE, &
                     t9ppg, "temperature", SE_DOUBLE, dppg, "dcoeff", SE_DOUBLE, &
                     yps_restart, "iso_massf", SE_DOUBLE_2D, m, nvar)
               call FSE_WRITE_IATTR(FID2, icountmodw, "shellnb",     m           )
               call FSE_WRITE_DATTR(FID2, icountmodw, "age",         age         )
               call FSE_WRITE_DATTR(FID2, icountmodw, "deltat",      dzeit       )
               call FSE_WRITE_DATTR(FID2, icountmodw, "total_mass",  xmrmax      )
               call FSE_WRITE_IATTR(FID2, icountmodw, "cycle",       icountmodw  )
               call FSE_CLOSE(FID2)

               open(restart_log, file='last_restart.out', status='replace')
               ! CR: for sub-time step output, write number of time (incl. sub-time) steps
               if ( ( isubmax  >  1 ) .and. ( writeoutsubcycle  ==  1 ) ) then
                  write(restart_log,*) icountmodw, icountmod, ipacket
               else
                  write(restart_log,*) icountmod, ipacket
               end if
               close(restart_log)
               deallocate (yps_restart)
            end if

            1337 if (iolevel >= 4) print *, "done output USEEPP .."

            ! write summary information
            if (iolevel >= 1) then
               idatastepsize = 10
               do j = 1, maxzonselt, idatastepsize
                  write(summary_output,9832) j, aITER(j), atiters(j), ansubt(j), antime(j), &
                        anuctime(j),anvar1(j) ,artime(j)
               end do
               write(summary_output,'(2(A21,ES10.3))') "# sum(solver_time) = ", &
                     sum(anuctime(1:maxzonselt)), "sum(physics_time) = ", &
                     sum(artime(1:maxzonselt))
               write(summary_output,'(A18,ES10.3)') "# mixing time/s = ", mixing_time
            end if
            9832 format(I4,X,I6,X,I6,X,I6,X,ES15.5,1X,ES13.3,X,I12,X,ES13.5)

            ! Check to see if the program has operated past the tstop parameter.
            ! If so, the program will close.  However, there are slaves, so the
            ! master needs to tell the slaves to finalize.  This is done by broadcasting
            ! the terminate_flag variable, which is .true. if the program is to terminate.
            ! If we are not stopping, all processors will just reach the 'end do' statements
            ! and execute the next iteration of the program

            call cpu_time(tcurrent)
            ttdiff = (tcurrent - tinit) / 3600._r8
            if (ttdiff > tstop) then ! Program terminates
               terminate_flag = .true.
            endif

            if (iolevel >= 4) then
               !     print some debugging output for the tstop parameter information
               !     every timestep         
               write(*,'(A9,i5,A17,L)') "Timestep:", icountmod, "terminate flag?", terminate_flag
               write(*,'(A6,ES12.3,A12,ES12.3,A11,F5.3)') "Time/s:", &
                     tcurrent,"tdiff/hrs= ", ttdiff, "tstop/hrs:", tstop
            end if
         end if
         ! end if (master)

         !     If terminate_flag = .True. then all processors terminate.
         !     Otherwise, they continue onwards to the next timestep or sub-timestep
         call broadcast(terminate_flag)

         if (terminate_flag) goto 1216

      end do subcycle
   end do cycleloop


   1216 if (master) then
      call cpu_time(tfin)
      if (iolevel  >=  1) then
         if (terminate_flag) then
            write(*,*) "Program terminated due to time limit tstop reached."
         else
            write(*,*) "Post-processing complete."
            write(summary_output,*) "Post-processing complete."
         end if
         write(*,*) "Total program time is:", tfin - tinit, "seconds"
         write(*,*) "Number of timesteps calculated:", icountmod - modstart
         write(summary_output,*)"Total program time is:",tfin - tinit, "seconds"
         write(summary_output,*) "Number of timesteps calculated:", icountmod - modstart
      end if
   end if

   close(summary_output)

   call mpi_finalize(ierr)


contains


   subroutine readhifinput(ing,ingperiod,xingest)
         ! read ppn_hif.input
         implicit none
         integer :: ing, ingperiod
         real(r8) :: xingest
         NAMELIST/ppn_hif/ ing,ingperiod,xingest

         ! defaults (from Falks RUN48)
         ing       = 1 ; ingperiod = 6 ; xingest   = 5.e-4_r8

         if ( master ) then
            open  ( 334 , file = 'ppn_hif.input' ) 
            read  ( 334 , NML  = ppn_hif         ) 
            close ( 334 ) 
         end if
   end subroutine readhifinput


   !> read ppn_frame.input
   subroutine readframeinput( t9threshold, ythreshold, modstart, modstop, &
               igrid, dxm, xmrmin, xmrmax, datdir, prefix, &
               nprno, nprnr, ioutc, trestart, tstop, &
               hif, sig_term_limit, isubmax, writeoutsubcycle )

         use communication
         use array_sizes
         use frame_knobs
         use utils

         implicit none

         integer :: fh, modstart, modstop, igrid, nprno, nprnr, ioutc, hif, isubmax, writeoutsubcycle
         real(r8) :: t9threshold, ythreshold, dxm, xmrmin, xmrmax, trestart, tstop, sig_term_limit
         character(len=256) :: datdir,prefix,ini_file_name2
         NAMELIST/ppn_frame/ iabuini, ini_filename, iolevel, &
               t9threshold, ythreshold, modstart, modstop, igrid, dxm, xmrmin, &
               xmrmax, code_source, datdir, prefix, nprno, nprnr, ioutc, &
               trestart, tstop, hif,  sig_term_limit, &
               isubmax, writeoutsubcycle, ini_filename, ini_filename2, ye_initial

         iabuini = 4            ! initialisation
         ini_filename = 'ini_abund.dat' ! if iabuini = 10 then use this file for 
         ! the initial abundance, if iabuini = 20 this file is 
         ! the abundance of the accreted material
         ini_filename2 = '../USEEPP/wd_abund.ppn' ! if iabuini = 20 file name for 
         ! the abundance distribution of the white dwarf
         iolevel = 3            ! how much output do you want, >4 is for
         ! debugging

         t9threshold = 0.008_r8
         ythreshold  = 1.d-99

         ye_initial = 0.5d0

         modstart = 1           ! start model for post-processing
         modstop  = 100

         igrid  = 2             ! set grid option
         dxm    = 1.d-2         ! grid size for static grid

         xmrmin = 0._r8            ! min mass coordinate for pp
         xmrmax = 2.0_r8           ! max mass coordinate

         code_source = 'MES'    ! GNV - EVL
         datdir = '/rpod2/fherwig/SEE/Set1.1/Set1.1-E0146'
         prefix = 'M2.00Z0.010'

         trestart = 1.5_r8         ! write restart file every trestart hours
         tstop    = 5000.       ! stop the program after tstop hours
         nprno  = 20            ! cycle interval for long output
         nprnr  = 50            ! cycle interval for restart outpu
         ioutc  = 500           ! number of cycles in output h5 files

         hif = 0                ! HIF implementation switch on(1)/off(0)

         sig_term_limit = 1d+10

         isubmax = 0            !switch on(1)/off(0) sub-time step implementation
         writeoutsubcycle = 0   !if 1, write out sub-time steps in hdf5 files (needs isubmax = 1)

         if ( master ) then
            fh = get_unused_file_handle()
            open ( fh , file = 'ppn_frame.input' ) 
            read ( fh , NML  = ppn_frame         ) 
            close( fh )

            if ( igrid == 1 ) then
               print *,'equidistand grid need to be done correctly'
               print *,'work in progress - Falk and Marco'
               stop
            end if    
         end if

         call broadcast(modstart)
         call broadcast(iolevel)
         call broadcast(isubmax)
         call broadcast(hif)

         isubmax = max(1, isubmax + 1)

   end subroutine readframeinput



   subroutine zoneselect(t9threshold,t9ppg,m,zone_to_compute,mmax)
         use array_sizes
         use utils, only: r8

         implicit none

         integer :: m, j, mmax
         logical :: zone_to_compute(msl)
         real(r8) :: t9ppg(msl), t9threshold
         mmax = 0      
         zone_to_compute(:) = .FALSE.
         do j = 1, m
            if (t9ppg(j) > t9threshold) then
               zone_to_compute(j)=.true.
               mmax = j
            end if
         end do
   end subroutine zoneselect



   subroutine master_pack_and_send( rho, t9, y0, xbuffer, buffer_size, slave_id, tag )
         use communication
         use utils, only: r8
         use array_sizes, only: nsp

         integer :: slave_id, tag, iposition, ierr
         integer(kind=8) :: buffer_size
         real(r8) :: rho, t9, y0(nsp), xbuffer(buffer_size)


         iposition = 0
         call mpi_pack(rho,       1, mpi_double_precision, xbuffer, buffer_size, iposition, mpi_comm_world, ierr)
         call mpi_pack(t9,        1, mpi_double_precision, xbuffer, buffer_size, iposition, mpi_comm_world, ierr)
         call mpi_pack(y0, size(y0), mpi_double_precision, xbuffer, buffer_size, iposition, mpi_comm_world, ierr)

         call mpi_send(xbuffer, buffer_size, mpi_packed, slave_id, tag, mpi_comm_world, ierr)
   end subroutine master_pack_and_send



   subroutine slave_pack_and_send( y0, xbuffer2, buffer_size2, tag )
         use communication
         use utils, only: r8
         use array_sizes, only: nsp

         integer :: iter, titers, nsubt, nvar1, tag
         integer(kind=8) :: buffer_size2
         real(r8) :: xbuffer2(buffer_size2), ntime, rtime, nuctime, y0(nsp)

         iposition = 0
         call mpi_pack(y0,     nsp, mpi_double_precision, xbuffer2, buffer_size2, iposition, mpi_comm_world, ierr)

         call mpi_send(xbuffer2, buffer_size2, mpi_packed, 0, tag, mpi_comm_world, ierr)
   end subroutine slave_pack_and_send



   subroutine master_recv_and_unpack(xbuffer2, buffer_size2, stat, y0)
         use communication
         use utils, only: r8
         use array_sizes, only: nsp

         integer :: stat(mpi_status_size), iposition, ierr
         integer(kind=8) :: buffer_size2
         real(r8) :: xbuffer2(buffer_size2), y0(nsp)

         iposition = 0

         call mpi_recv(xbuffer2, buffer_size2, mpi_packed, mpi_any_source, mpi_any_tag, mpi_comm_world, stat, ierr)

         call mpi_unpack(xbuffer2, buffer_size2, iposition, y0,      size(y0),  mpi_double_precision,     mpi_comm_world, ierr)
   end subroutine master_recv_and_unpack


   subroutine slave_recv_and_unpack(rho, t9, y0, xbuffer, buffer_size, tag)
         use communication
         use array_sizes, only: nsp
         use utils, only: r8

         integer(kind=8) :: buffer_size
         real(r8) :: rho, t9, y0(nsp), xbuffer(buffer_size)
         integer :: stat(mpi_status_size), iposition, ierr, tag

         call mpi_recv(xbuffer, buffer_size, mpi_packed, 0, mpi_any_tag, mpi_comm_world, stat, ierr)
         tag = stat(mpi_tag)
         if (tag == 0) then
            return
         end if
         iposition = 0
         call mpi_unpack(xbuffer, buffer_size, iposition, rho,   1, mpi_double_precision, mpi_comm_world, ierr)
         call mpi_unpack(xbuffer, buffer_size, iposition,  t9,   1, mpi_double_precision, mpi_comm_world, ierr)
         call mpi_unpack(xbuffer, buffer_size, iposition,  y0, nsp, mpi_double_precision, mpi_comm_world, ierr)
   end subroutine slave_recv_and_unpack



   subroutine terminate_slave(slave_id)
         use communication
         integer :: slave_id, ierr

         call mpi_send(1.0_r8, 0, mpi_integer, slave_id, 0, mpi_comm_world, ierr)
   end subroutine terminate_slave



   subroutine terminate_all_slaves()
         use communication
         integer :: i
         do i = 1, num_procs - 1
            call terminate_slave(i)
         end do
   end subroutine terminate_all_slaves



   subroutine iniabund( m, xm,y_accrete, mwd, yps, irestart, fileout, nvar )
         ! ^_^ see physics/abundances.F90 for definition of iniabund options
         use array_sizes
         use frame_knobs
         use nuc_data
         use utils
         use constants
         use abundances
         implicit none

         integer :: i, m, nvar, irestart, nvar_restart, ii, isomeric_state_restart(nsp), &
               izn_restart, ian_restart, iisomer, mwd, fid2

         character(len=256) :: fileout
         real(r8) :: yps(msl,nsp), y0(nsp), xm(msl), yps_restart(msl,nsp), y_accrete(nsp), &
               y0wd(nsp), an_restart(nsp), zn_restart(nsp), an, zn

         common / cnetw / an(nsp), zn(nsp)

         select case(iabuini)
         case default
            ! initialise uniform composition
            call set_initial_abundances(y0)
            do i = 1, m
               yps(i,:) = y0(:)
            end do

         case(20)
            ! ^_^ accreting WD model

            iabuini = 11
            ! accreted composition:
            ini_filename = ini_filename
            call set_initial_abundances(y0)
            do i = mwd + 1, m
               yps(i,:) = y0(:)
            end do

            ! WD abundance:
            ini_filename = ini_filename2
            call set_initial_abundances(y0wd)
            do i = 1, mwd
               yps(i,:) = y0wd(:)
            end do

            iabuini = 20

         case(0)
            ! ^_^ mppnp restart file
            write(*,*) "opening restart file ", fileout
            call check_file_exists(fileout)

            call FSE_OPEN(fileout,FID2)
            call FSE_READ_IATTR(FID2,-1,"zisnb", nvar_restart)

            if (nvar > nvar_restart) write(*,*) "WARNING: mppnp/iniabun: new species will be initialised to zero."

            call FSE_READ_IATTR(FID2, irestart, "shellnb", m)
            if (m > msl) then
               print *, "INIABUND: parameter MSL is too small.  M = ", m, "MSL = ", msl
               stop
            end if

            call FSE_READ_D(FID2, irestart, m, "mass", xm)
            call FSE_READ_D_2D(FID2, irestart, m, "iso_massf", yps_restart, msl,nvar_restart)
            call FSE_READ_DARRAYATTR ( FID2 , -1 , "Z"              , zn_restart             , nvar_restart ) 
            call FSE_READ_DARRAYATTR ( FID2 , -1 , "A"              , an_restart             , nvar_restart ) 
            call FSE_READ_IARRAYATTR ( FID2 , -1 , "isomeric_state" , isomeric_state_restart , nvar_restart ) 
            call FSE_CLOSE(FID2)

            !  WARNING: restarts should be done using same network, either by:
            !     1. using ininet = 3
            !        or
            !     2. same ppn_physics.input settings

            yps(:,:) = ZERO

            do i = 1, nvar_restart
               ian_restart = nint ( an_restart ( i ) ) 
               izn_restart = nint ( zn_restart ( i ) ) 
               iisomer     = isomeric_state_restart(i)
               ii          = niso(ian_restart,izn_restart,iisomer)
               if (ii == 0) then
                  write(*,*) "WARNING mppnp/iniabund: in restart requesting", &
                        " species that is not available in networksetup.txt."
                  write(*,*) ian_restart,izn_restart,iisomer
               else
                  if (nint(an(ii)).ne.ian_restart .or. nint(zn(ii)).ne.izn_restart .or. isomeric_state(ii).ne.iisomer) then
                     write(*,*) "mppnp/iabunini: a,z,isomer in restart can not be matched to database"
                     write(*,*)"a db vs restart:", nint(an(ii)), ian_restart 
                     write(*,*)"z db vs restart:", nint(zn(ii)), izn_restart 
                     write(*,*)"iso db vs restart:", isomeric_state(ii), iisomer
                     write(*,*)"BAD WARNING!!! only ok for NEUT"
                  endif
                  yps(1:m,ii)  = yps_restart(1:m,i)
               end if
            end do

            do i = 1, m
               if (abs(sum(yps(i,:)) - 1._r8) > SUM_TOL) then
                  write(*,*) "WARNING: zone: ",i, "sum(yps)=",sum(yps(i,:)), ' , renormalising '
                  yps(i,:) = yps(i,:)/sum(yps(i,:))
               end if
            end do 
         end select

         write(*,*) "finished initialising abundances"

         ! save composition of accreting material for further accretion
         y_accrete(:) = y0(:)
   end subroutine iniabund



   subroutine customgrid(xmrmin,xmrmax,dxm,xm,m)

         ! *** implement a simple static grid

         use frame_knobs
         use utils, only: r8

         implicit none

         integer :: m,i
         real(r8) :: xm(m), xmrmin, xmrmax, dxm

         ! *** in mppnp Lagrangian grid points are meant to be the location of 
         !     the center of the mass cell, which implies that the cell faces 
         !     are both dxm_i/2 away from xm(i)

         i = 1
         xm(1)=xmrmin+dxm/2.d0
         do while (xm(i)+dxm/2.d0 < xmrmax)
            i = i+1
            xm(i) = xm(i-1) + dxm
         end do
         m = i - 1

         if (iolevel >= 4) then
            print *,"Customgrid:"
            do i=1,m
               print *,i,xm(i)
            end do
         endif
   end subroutine customgrid



   subroutine customgrid_hif( xmrmin, xmrmax, dxm, xm, m )

         ! *** In this special version for post-processing the Sakurai evolution
         !     model we design our own statically refined grid according to
         !     BFI-HIF081210-1
         ! *** in mppnp Lagrangian grid points are meant to be the location of 
         !     the center of the mass cell, which implies that the cell faces 
         !     are both dxm_i/2 away from xm(i)
         ! *** WARNING: this subroutine is not debugged, the code below was just moved
         !     here in the major reshuffle of mppnp in July 2009. it may still work

         use frame_knobs

         implicit none

         integer m,i
         double precision xm(m),xmrmin,xmrmax,dxm


         XM(1) = xmrmin
         i = 1
         if (iolevel  >=  1) write(*,*) "xm(1), xmrmin",xm(1), xmrmin
         if (iolevel  >=  3) open(87,file="grid.DAT")


         BGRID:do
            if (xm(i)+5.d-4 < 0.60D0) then
               dxm=5.d-4
            else if (xm(i)+5.d-5 < 0.60417D0) then
               dxm=5.d-5
            else if (xm(i)+1.d-5 < 0.60427D0) then
               dxm=1.d-5
            else 
               dxm=1.d-6
            end if
            i = i+1
            xm(i) = xm(i-1) + dxm
            if ( iolevel  >=  3 ) write(87,*) i-1,xm(i-1)
            if (xm(i)  >=  xmrmax) then
               xm(i) = xmrmax
               exit BGRID
            end if
         end do BGRID
         M = i

         if ( iolevel  >=  1 ) write(*,*) "End of Grid setup: no grids = ", M, xm(M)
         if ( iolevel  >=  3 ) close(87)
   end subroutine customgrid_hif



   subroutine surf_elem( isurf, isoz, isoa, yps_surf, nvar, num_of_elements, &
               elem_surf, num_surf, elem_index, elem_num, dimens )

         ! *** This subroutine calculates the elemental abundances and returns the result
         !     as an array that can be printed into the surf.h5 file

         use array_sizes
         use utils, only: r8

         implicit none

         integer :: dimens  ! number of isotopes at any one value of Z
         real(r8) ::  isoz(nsp), isoa(nsp), yps_surf(nvar)
         real(r8) ::  elem_surf(num_of_elements)
         real(r8) ::  num_surf(num_of_elements)
         real(r8) ::  xsum1, xsum2
         integer :: i, j, k, isurf, num_of_elements, nvar
         integer :: elem_index(num_of_elements, dimens)
         integer :: elem_num(num_of_elements)
         !     elem_index contains the array indices for yps that contain elements of a
         !     particular Z.  For example, if elem_index(1,:) = (/2, 3, 7/) then the
         !     Hydrogen isotopes (Z=1) are contained in yps(:,2), yps(:,3) and yps(:,7).
         !     The summation will take the surface mass fraction abundances, so the final sum is:
         !     H = yps(m,2) + yps(m,3) + yps(m,7)

         !     elem_num contains the number of isotopes for each element.  This is to save
         !     some time so that the loop does not have to go over the whole of dimens
         !     and may be useful to know for debugging.

         ! *** The first step is to determine which isotopes are in the network and build
         !     up an index array containing isotopes of common z.  This is stored so the
         !     same step does not need to be repeated at every output (the network does
         !     not modify isotopes on the fly anyhow).

         !     isurf = 0:  calculate the index array for isotopes of a certain element
         !     isurf = 1:  already have it, so no need to do it again

         if (isurf  ==  0) then
            elem_index(:,:) = 0

            !     loop: neutrons, hydrogen, helium, etc
            do i = 0, num_of_elements-1
               k = 1
               do j = 1, nvar
                  if (nint(isoz(j))  ==  i) then
                     elem_index(i+1,k) = j
                     k = k+1
                  end if
               end do
               elem_num(i+1) = k - 1
            end do
            isurf = 1
         end if
         ! *** Now that the index array is generated, we now know which isotopes to sum up
         !     to determine the elemental abundances
         do i = 0, num_of_elements-1
            xsum1 = 0.d0
            xsum2 = 0.d0
            if (elem_num(i+1) .ne. 0) then
               do j = 1, elem_num(i+1)
                  k = elem_index(i+1,j)
                  xsum1 = xsum1 + yps_surf(k)
                  if (i .ne. 0) then
                     xsum2 = xsum2 + (yps_surf(k)/isoa(k))
                  end if
               end do
            end if
            elem_surf(i+1) = xsum1
            num_surf(i+1) = xsum2
         end do
   end subroutine surf_elem



   subroutine surf_aver(nshells1, nvar, xm, yps, surf_coord, yps_av)
         use frame_knobs
         use utils, only: r8

         implicit none

         real(r8) :: xm(nshells1)
         real(r8) :: xm_temp1(nshells1), xm_temp2(nshells1+1)
         real(r8) :: yps(nshells1, nvar), yps_int(nvar), yps_av(nvar)
         real(r8) :: yps_temp1(nshells1, nvar)
         real(r8) :: yps_temp2(nshells1+1, nvar)
         real(r8) :: surf_coord

         real(r8) :: DeltaX(nvar), Deltaxm, xarea
         integer :: i, ishell, nshells1, nshells2, nvar
         !     nshells1: total number of shells in xm, yps, etc
         !     nshells2: total number of mesh points to integrate over
         !     (calculated later)

         ! *** This subroutine performs an average of the abundances at the
         ! surface according
         !     to the input parameters surf_coord

         ! *** The subroutine operates on the pretense that the mass coordinate
         ! increases
         !     with shell number.  If this is not the case, it reverses the
         !     direction
         if ( xm(1)  >  xm(2) ) then
            do i = 1, nshells1
               xm_temp1(i) = xm(nshells1+1-i)
            end do
            xm(:) = xm_temp1(:)
         end if
         ! *** There will be some coordinate which defines where to average to.
         ! Therefore,
         !     an interpolation is required to determine the abundances at that
         !     point.  This
         !     is just a simple straight line interpolation in log space.
         if (surf_coord  <=  xm(1)) then
            !        surf_coord is at the core (average over the whole star)
            if (iolevel  >=  3) then
               print *,"Surf_elem: surf_coord is at the core."
               print *,"Averaging over the whole star..."
            end if
            yps_temp1(:,:) = log10(yps(:,:))

            Deltaxm = xm(nshells1) - xm(1)
            do i = 1, nvar
               call trapz(xm, yps_temp1(:,i), nshells1, xarea)
               yps_av(i) = xarea / Deltaxm
            end do
            return

         else if (surf_coord  >=  xm(nshells1)) then
            ! surf_coord is at the surface; just take the surface abundances
            if (iolevel  >=  3) then
               print *,"Surf_elem: surf_coord is at or beyond the surface."
               print *,"Taking surface abundances..."
            end if
            yps_av(:) = yps(nshells1,:)
            return
         else
            do i = 1, nshells1-1
               if (xm(i) < surf_coord .and. xm(i+1) >= surf_coord) then
                  ! surf_coord is between xm(i) and xm(i+1): interpolate
                  DeltaX(:) = (log10(yps(i+1,:)) - log10(yps(i,:))) / (xm(i+1)-xm(i))
                  yps_int(:) = log10(yps(i,:)) + DeltaX*(surf_coord-xm(i))
                  ishell = i + 1
                  exit
               end if
            end do
            xm_temp2(1) = surf_coord
            yps_temp2(1,:) = yps_int(:)
            if (ishell  ==  nshells1) then
               xm_temp2(2) = xm(ishell)
               yps_temp2(2,:) = log10(yps(ishell,:))
            else
               xm_temp2(2:) = xm(ishell:)
               yps_temp2(2:,:) = log10(yps(ishell:,:))
            end if
            nshells2 = nshells1 - ishell + 2
            Deltaxm = xm_temp2(nshells2) - xm_temp2(1)
            do i = 1, nvar
               call trapz(xm_temp2, yps_temp2(:,i), nshells2, xarea)
               yps_av(i) = xarea / Deltaxm
            end do
            return
         end if
   end subroutine surf_aver



   subroutine trapz(x, y, n, xarea)
         implicit none

         integer i, n
         double precision x(n), y(n), xarea

         ! *** This routine performs trapezoidal integration: it sums up areas of trapezes
         !     over the whole x range and outputs the area under the curve
         xarea = 0.d0
         do i = 1, n-1
            xarea = xarea + (0.5_r8*(y(i)+y(i+1))*(x(i+1)-x(i)))
         end do
   end subroutine trapz



   subroutine decay_I( t9, nvar, age, nvrel, v_terr, yps_surf_nodecay, yps_surf_decay )

         ! *** This routine performs decay of the yps array in the argument
         ! ***  (yps_surf_nodecay), providing new array that considers the decay
         ! *** contribution (yps_surf_decay).

         use array_sizes
         use nuc_data
         use backward_euler
         use physics_knobs
         use frame_knobs
         use rates
         use utils, only: r8, calculate_ye

         implicit none

         integer :: i, nvar, nvrel, ierr
         real(r8) :: t9, dzeit_decay, rho_fake, age, ye, an(nsp), zn(nsp), &
               yps_surf_nodecay(nsp), yps_surf_decay(nsp), v_terr(nre), t9_fake

         common/cnetw/an,zn

         rho_fake = 1._r8
         call evaluate_all_rates( ye, nvar, nvrel, rho_fake, t9, yps_surf_nodecay, nu )

         ! *** dzeit_decay = time step to have complete decay.
         dzeit_decay = 3.1536d+16  ! 1Gyr 

         call calculate_ye(yps_surf_nodecay, an, zn, ye, considerisotope)

         yps_surf_decay(:)  = yps_surf_nodecay(:)

         !> this probably isn't a good way, or at least needs rethinking
         call integrate_network( nvar, yps_surf_decay, t9, t9, rho_fake, rho_fake, ye, dzeit_decay, nvrel, nu, ierr )

         ! *** at this point I have calculated the decayed abundances.
   end subroutine decay_I



   subroutine opening_statement()
         write(*,*) " This is MPPNP: multi-zone post-processing network."
         write(*,*) " (c) 2008 Marco Pignatari, Falk Herwig,"
         write(*,*) "     Michael Bennett, Raphael Hirschi and the"     
         write(*,*) "     NuGrid collaboration."
         write(*,*)
   end subroutine opening_statement

end program mppnp


function ylin(x,x1,x2,y1,y2)
      use utils, only: r8
      real(r8) :: ylin, x, x1, x2, y1, y2

      ylin=y1+((y2-y1)/(x2-x1))*(x-x1)
end function ylin







! *** variable definitions *******************************************************
!    nsp           max number of isotopes 
!    msl           max number of spatial grid zones for multi-zone 
!                  calculations (not supported in this version)
!    nre           max number of reactions that can be read in, this parameter 
!                  should be larger than nnn+nvcp in parameter.input
!    nvrel         number of used reactions
!    v             reaction rate array
!    yps           the master abundance array [mass fractions]
!    zis           character array that holds the look-up table of isotope 
!                  names
!                  function ispe returns index number of zis array for 
!                  species name
!    considerisotope(nsp),considerreaction(nre)
!                  these two logical arrays control if an isotope or a 
!                  reaction is actually included in the network calculation
!    dzeitj        time step in years (zeit=time, jahre=years)
!    dzeit         time step in seconds
!    iabuini       how are initial abundances set, see subroutine iniabund
!    modstart      start model
!    modstop       stop model
!    t9threshold   activate network only for T9 above this threshold
!    ythreshold    if abundance lower than ythresh then abundance = ythresh
!    zn            proton numbers of isotopes used in yps
!    an            mass numbers of isotopes used in yps
!    isomeric_state a 1D integer array containing isomeric states of isotopes
!    mwd           mesh point at which initial abundances change from the envelope to wd-core ones
!    xmwd          mass at which initial abundances change from the envelope to wd-core ones
!    yinp          accretion composition(?)
!    y_accrete     accretion composition(?)
! *** 1D post-processing: *********************************************
!    maxmod        max number of models to be processed (just for assigning arrays)
!    maxfile       max number of files to be processed (just for assigning arrays)
!    age           age in stellar evolution sequence, only read and written, in yrs
!    xmcoord       mass coordinates of 1: surface
!                  in stellar          2: H-shell
!                  evolution model     3: He-shell
!    code_source  which code was used for the stellar structure?'GNV','MES','EVL'
!    datdir       character name of path to DATEN directory with HPF files
!    prefix       prefix for the input hdf5 file (and the index .idx file as well)
!    filein       full name for the input hdf5 file
!    fileout      full name for the output hdf5 file
!    msl          max number of cells in pp-grid, >32000 => needs -i8
!    xmm(msl)    mass grid in input SE
!    dse(msl)    D (diffusion coefficient) from stellar evolution (SE)
!    rhose(msl)  rho from stellar evolution (SE)
!    t9se(msl)   T9 from stellar evolution (SE)
!    rse(msl)    R from SE, in Rsun
!    xm(msl)      pp-grid, Lagrangian, in solar masses
!    dppg(msl)    D (diffusion coefficient) interpolated on pp-grid
!    sig_term_limit a limit for diffusion coefficients similar to that in mesa
!    rhoppg(msl)  rho on pp-grid
!    t9ppg(msl)   T9 on pp-grid
!    rppg(msl)    R on pp-grid
!    xmrmax       max mass of pp-grid
!    xmrmin       min mass of pp-grid
!    dxm          delta_m in amr base grid (iamr(1,:)=1
!    m            number of grid points of pp-grid, determined during grid setup
!    icountmod    the cycle number in the stellar evolution input and in mppnp
!    dconst       4*Pi*Rsun^2/Msun for calculating D_Langrangian from D_Eularian
!    nprnr        cycle interval for restart cycles
!    nprno        cycle interval for standard output (ABU files, depricated)
!    aITER        Summary information: number of iterations for each shell
!    antime       Summary information: time taken f<or NRNW function call
!    ansubt       Summary information: number of timesteps in nucnet99 (incl. subtime)
!    artime       Summary information: time taken for rnetw2008 function call
!                 (or testnse)
!    anuctime     Summary information: time taken for nucnet99 function call
!    idatastepsize This variable determines how much output information is generated
!    isurf        Switch that determines whether to determine the isotope index array
!                 for elements.  It is used to calculate elemental yields.
!    elem_surf    array containing element abundances at the surface
!    xnum_surf    array containing number abundances of elements at the surface
!    num_of_elements number of elements being considered in network (including neutrons)
!
! *****************************************************************

! *********************************************************************
   !     nnn                    : Number of isotopes in the network.
   !     i_nse                  : Flag for setting up nse functionality.
   !     ilastcycle             : Slaves follow the cycleloop and 
   !                              therefore need to know the maximum 
   !                              number of models being calculated.
   !     istart                 : Tables are read when istart=0, therefore
   !                              slaves need to be told that istart=1.
   !     nvar                   : Number of considered rates as defined in        
   !                              rnetw2007.Used by slaves in rnetw2008.
   !     nvrel                  : Number of considered reactions as defined
   !                              in rnetw2007.Used by slaves in rnetw2008.
   !     ininet                 : Input parameter to determine the network
   !                              setup in rnetw2007 and rnetw2008.
   !     index_reaclib          : Integer that determines which version of
   !                              reaclib to use: basel or JINA
   !     ztot                   : NetworkI input data: proton numbers as
   !                              read from isotopedatabase.txt.
   !     atot                   : NetworkI input data: proton numbers as
   !     v                      : Used as an input argument for rnetw2008.
   !                              This is needed for the first function
   !                              call.
   !     an                     : ppn_physics.input: mass numbers.
   !     zn                     : ppn_physics.input: proton numbers.
   !     amin                   : NetworkI input data: minimum mass of an
   !                              isotope as read from isotopedatabase.txt
   !     eltot                  : NetworkI input data: Isotope array from
   !                              isotopedatabase.txt.
   ! *********************************************************************


   ! *** MPI declarations and function calls to set up the parallel environment**********
   !	buffer_size		   : size of the buffer of message from Master -> Slaves
   !  buffer_size2      : size of the buffer of message from Slaves -> Master
   !	xbuffer			   : message buffer Master -> Slaves
   !  xbuffer2          : message buffer Slaves -> Master
   !	num_sent		      : Total number of messages sent.
   !  irecv             : Number of messages to receive.
   !	iposition		   : Position of the pointer used in the function MPI_PACK.
   !	isender			   : Stores which process sent the message.
   !	tag			      : message tag used to store which meshpoint abundances
   !                                 are being calculated for. 
   !	imass_shell		   : This stores the meshpoint that the answer, y0, will 
   !                                 eventually be placed in yps (from tag).
   !  stat(MPI_STATUS_SIZE)   : A handy MPI variable that stores information on the status
   !             of a message such as the tag, source, ...
   ! ***********************************************************************************
