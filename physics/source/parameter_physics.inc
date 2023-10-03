! *** delclare parameters related to physics package

! *** in ppn_physics.input one may run with
! (FH@MP,141130: please check this - can we still run with NNN=0?)
! *** for NNN=0   - needs to be updated because of VITAL network extension
!      integer, parameter :: NSP=36, NRE=65 ! needs to reflect choices in 
                                            ! ppn_physics.input
! *** for NNN > 0 - running with additional background network                
!      integer, parameter :: NSP=1150, NRE=15000 ! isotopedatabase.txt
      integer, parameter :: NSP=5250, NRE=75000 ! isotopedatabase_all.txt
!      integer, parameter :: NSP=5640, NRE=78000 ! isotopedatabase_all.txt
!      integer, parameter :: NSP=300,  NRE=15000 ! Pavel's nova case

! *** choices for reaclib
! *** reaclib parameters --> Basel
! *** INDEX_REACLIB = 0
!      integer, parameter :: ndim_rl=61721, nnpartdim=5427
! *** reaclib parameters --> Jina (without species beyond At):
! *** INDEX_REACLIB = 1, 2
      integer, parameter :: ndim_rl=74313, nnpartdim=5497
! *** INDEX_REACLIB = 3 -- JINA raclib extended to Cf --
!      integer, parameter :: ndim_rl=77953, nnpartdim=6789
