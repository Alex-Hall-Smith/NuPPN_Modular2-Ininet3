module output
   use utils
   use nuc_data, only: zis
   use array_sizes, only: nsp
   implicit none

   integer, public :: xtime_fh
   character(len=80), public :: xtime_cycfmt, xtime_hedfmt
   character(len=10), parameter, public :: xtime_fn = 'x-time.dat'
   character(len=33), parameter, public :: stdout_fmt = "(I6,8(1x,1PD9.1),1X,I5,3(1x,I4))"

contains

   subroutine output_init(nvar,iprint_true)
         ! *** write header and initial condition
         integer :: i, num_cols, nvar, iprint_true(nsp)

         xtime_fh = get_unused_file_handle()
         num_cols = nvar + 4  ! for printing x-time.dat columns
         write( xtime_cycfmt, '("(I7,1x1es12.5,1x,f8.4,", I0, "(1x,1ES12.5))")' ) num_cols
         write( xtime_hedfmt, '("(5a,", I0 ,"(A1,I4,A1,a5,2x))")' ) num_cols

         open( xtime_fh, file = xtime_fn )

         write( xtime_fh, xtime_hedfmt ) '#|cycle |time         ', '|t9     ', '|rho         ', &
               '|1-sum(yps)' , '  |ye         ', ('|', i + 6, '-', zis(iprint_true(i)), i = 1, nvar)

         ! *** STDOUT
         write(*,'(3A)') ' cycle   age       N_n       T_9 ' // &
               '      rho       ye        <tNRNW>   tN_last' // &
               '   tminv_l Nspec   IT  TIT nsubt'

   end subroutine output_init



   subroutine output_fin()
         close(xtime_fh)
   end subroutine output_fin

end module output
