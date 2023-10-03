module jbj16
   use frame_knobs, only: iolevel
   use physics_knobs, only: jbj_mode
   use utils, only: r8, irate, idrdt, idrdd, get_unused_file_handle, generic_interp_2d
   use constants
#ifndef PPN
   use communication
#endif
   implicit none
   private

   ! ^_^ Weak reaction rates (beta-, beta+, e+-capture, e--capture) from Jones, Bertolli & Johnson (submitted 2016).
   !     To use these rates, you will need the raw data files that are not available from the git server that hosts the code, owing
   !     to their size. They will be available on the NuGrid/CANFAR CADC data server once published. For now, request access from
   !     Sam. If you do not have the data, you can still run the code as normal by setting jbj_mod = 1 in pn_physics.input (default)

   ! ^_^ CONTROLS:
   !     jbj mode is set in ppn_physics,input (see physics_knobs.F90 in the physics package)
   !     since JBJ16 rates cover the entire region of the isotopic chart that the Oda+ '94 rates cover, you can either:
   !        jbj_mode = 1 : don't use JBJ rates at all, just use Oda rates
   !        jbj_mode = 2 : use JBJ rates only where Oda rates are used, i.e. Oda rates are directly replaced
   !        jbj_mode = 3 : use full set of rates from Oda+ '94, and add remaining rates from JBJ where there is no overlap 
   !        jbj_mode = 4 : use JBJ rates everywhere provided, thus using none of the Oda rates (and extending beyond)

   character(len=256), parameter :: listofrates = '../NPDATA/jbj16/_ratelist.txt'
   integer, parameter :: &
         num_jbj_rates = 288, &
         nT = 145, & ! temperature (T9) dimension
         nr = 501, & ! electron density dimension
         jbjminA = 17, jbjmaxA = 40, jbjminZ = 8, jbjmaxZ = 20
   character(len=4) :: weaktype(num_jbj_rates) ! 'beta' or 'ecap'
   real(r8) :: &
         T9(nT), &                           ! temperature coordinates
         log10T(nT), &                       ! log10 temperature coordinates
         lyr(nr), &                          ! electron density coordinates
         raw_rates(num_jbj_rates,nr,nT), &   ! the whole grid of reaction rates
         jbj_rate(num_jbj_rates), &          ! holds the rate for each reaction(3) for the required T and rhoye
   dlT, dT9, dlyr, &                   ! grid spacings
         minT9, maxT9, minlyr, maxlyr, minlT, maxlT   ! grid bounds
   integer :: &
         Ai(num_jbj_rates), Zi(num_jbj_rates), &   ! parent nucleus A and Z for each reaction
         Af(num_jbj_rates), Zf(num_jbj_rates), &   ! daughter nucleus A and Z for each reaction
         jbj_ntrans( jbjminA:jbjmaxA, jbjminZ:jbjmaxZ, 1, 13:14 ) ! (A, Z, 'state', reaction type); for merging into main network
   ! reaction type 13 = b- or e+-capture; 14 = b+ or e--capture
   public jbj_init, jbj_rate, jbjminA, jbjmaxA, jbjminZ, jbjmaxZ, jbj_ntrans, jbj_interpolate_rates

contains



   subroutine jbj_init()
         character(len=256) :: rate_fn
         integer :: i

         if (jbj_mode == 1) return

         jbj_ntrans(:,:,:,:) = 0

#ifndef PPN
         if (ipid == 0) then
#endif
            call jbj_read_data()
#ifndef PPN
         end if
         call jbj_broadcasts
#endif

         log10T = log10( T9 * 1.e9_r8 )
         do i = 1, num_jbj_rates
            if ( Zf(i) > Zi(i) ) then
               weaktype(i) = 'beta'
               jbj_ntrans( Ai(i), Zi(i), 1, 13 ) = i
            else
               weaktype(i) = 'ecap'
               jbj_ntrans( Ai(i), Zi(i), 1, 14 ) = i
            end if
         end do

         ! ^_^ tables are uniformly spaced in log(T / K) and log(rho * Ye)
         dlT   = sum(log10T(2:) - log10T(:nT-1)) / (nT - 1._r8)
         dlyr  = sum(lyr   (2:) - lyr   (:nr-1)) / (nr - 1._r8)

         minT9    = minval(T9);  maxT9    = maxval(T9)
         minlyr   = minval(lyr); maxlyr   = maxval(lyr)
         minlT    = log10( minT9 * 1.e9_r8 ); maxlT   = log10( maxT9 * 1.e9_r8 )

   end subroutine jbj_init



   subroutine jbj_interpolate_rates(ye_in, rho_in, T9_in)
         ! ^_^ perform bilinear interpolation (in log10 electron density and temperature in GK) of reaction rates
         real(r8), intent(in) :: ye_in, rho_in, T9_in ! electron density and temperature (GK) at which the rates are needed
         real(r8) :: lyr_in, myT9, mylyr, mylT, &
               xa, xb, ya, yb, &       ! interpolation variables
               T91, T92, lyr1, lyr2, & ! bounding values of the interpolation
               f11, f21, f12, f22, &   ! bounding data values
               fA, fB, fC, fD
         real(r8), dimension(num_jbj_rates) :: lr, dlrdt9, dlrdlyr
         integer :: i, &
               idxlT, idxT9, idxlyr           ! indices of lower bounding values

         lyr_in = log10( ye_in * rho_in )

         ! clip values to table edge
         myT9  = min( maxT9,  max( T9_in,  minT9  ) )
         mylyr = min( maxlyr, max( lyr_in, minlyr ) )
         mylT  = log10( myT9 * 1.e9_r8 )

         ! ^_^ get index using grid spacing (constant spacing in log(T / K) and log(rho * Ye) )
         idxlT    = floor( ( mylT  - minlT  ) / dlT  ) + 1
         idxT9    = idxlT
         idxlyr   = floor( ( mylyr - minlyr ) / dlyr ) + 1

         call generic_interp_2d(nr, nT, lyr, T9, num_jbj_rates, raw_rates(:,:,:), mylyr, myT9, &
               idxlyr, idxlyr+1, idxT9, idxT9+1, lr, dlrdt9, dlrdlyr)

         jbj_rate(:) = 10._r8 ** lr

         if (.false.) then
            print *, '************************************'
            print *, 'you asked for rate at T9, logYeRho =',T9_in, lyr_in, myT9, mylyr, mylT
            print *, 'this was found to be between temperatures', T9(idxT9), T9(idxT9+1) 
            print *, 'and densities', lyr(idxlyr), lyr(idxlyr+1)
            print *, 'with indices', idxT9, idxT9+1, idxlyr, idxlyr+1
            print *, 'bounding values of last rate:', f11, f21, f12, f22
            print *, 'interpolation constants:', xa, xb, ya, yb, fA, fB, fC, fD
            print *, 'table limits:', minT9, maxT9, minlyr, maxlyr
            print *, '************************************'
         end if

   end subroutine jbj_interpolate_rates



   subroutine jbj_read_data()
         character(len=256) :: rate_fn
         integer :: i, ierr, Tdim, rdim, rate_list, rate_data
         real(r8) :: temp_array(nr, nT)

         if ( iolevel >= 1) print *, 'reading weak rates from JBJ16...'
         rate_list = get_unused_file_handle()
         open( rate_list, file = listofrates, action = 'read' )
         do i = 1, num_jbj_rates
            ! ^_^ read file name containing rate data
            read( rate_list, *, iostat = ierr ) rate_fn
            if ( ierr /= 0 ) then
               print *, 'error reading file:', trim(listofrates), '; iostat = ', ierr
               stop
            end if
            rate_fn = '../NPDATA/jbj16/' // trim( rate_fn )
            ! read rate data from that file
            if ( i == 1 ) rate_data = get_unused_file_handle()
            open( rate_data, file = rate_fn, form = 'unformatted', access='stream', action = 'read' )
            read( rate_data ) Ai(i), Zi(i)
            read( rate_data ) Af(i), Zf(i)
            read( rate_data ) rdim, Tdim
            if (Tdim /= nT .or. rdim /= nr) then
               print *, 'unexpected table dimensions for jbj16 rate:', rate_fn
               close( rate_data )
               close( rate_list )
               stop
            end if
            read( rate_data ) lyr(:); read( rate_data ) T9(:)
            ! ^_^ read primary rate
            read( rate_data ) raw_rates(i,:,:)
            read( rate_data ) temp_array
            ! ^_^ read complimentary rate and add to total rate
            raw_rates(i,:,:) = log10( 10._r8 ** raw_rates(i,:,:) + 10._r8 ** temp_array(:,:) )
            close( rate_data )
         end do

         close( rate_list )

   end subroutine jbj_read_data



#ifndef PPN
   subroutine jbj_broadcasts
         call broadcast(raw_rates)
         call broadcast(Zi)
         call broadcast(Zf)
         call broadcast(Ai)
         call broadcast(Af)
         call broadcast(lyr)
         call broadcast(T9)
   end subroutine jbj_broadcasts
#endif


end module jbj16
