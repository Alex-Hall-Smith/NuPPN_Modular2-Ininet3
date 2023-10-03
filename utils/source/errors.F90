module errors
   implicit none
   private

   integer, parameter, public :: &
         SLV_CNVRG = 1, &
         MAX_SUBSTEPS = SLV_CNVRG + 1, &
         BAD_STEP = MAX_SUBSTEPS + 1, &
         NEG_ABUNDS = BAD_STEP + 1, &
         HUGE_ABUNDS = NEG_ABUNDS + 1, &
         SUM_X = HUGE_ABUNDS + 1, &
         GOT_NAN = SUM_X + 1, &
         BAD_YE = GOT_NAN + 1, &
         NSE_CNVRG = BAD_YE + 1, &
         NSE_WEAK = NSE_CNVRG + 1, &
         SMALL_DT = NSE_WEAK + 1

   character(len=11), dimension(11), parameter, public :: &
         error_code = [ 'CONVERGENCE', '  MAX_STEPS', '   BAD_STEP', &
         ' NEGATIVE_X', '     HUGE_X', '      SUM_X', '        NAN', '     BAD_YE', &
         '  NSE_CNVRG', '   NSE_WEAK', 'DT_UNDRFLOW' ]

end module errors
