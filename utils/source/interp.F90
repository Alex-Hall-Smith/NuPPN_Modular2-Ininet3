module interp
   use utils, only: r8, bsearch_r8
   implicit none
   private

   type, public :: interpolant
      integer :: siz
      real(r8), allocatable :: x(:), y(:)
   end type interpolant


   public create_1d_interpolant, interpolate_1d

contains

   !> create an interpolating function for a set of 1d data
   subroutine create_1d_interpolant( siz, x, y, i )
         integer :: siz
         real(r8), intent(in) :: x(siz), y(siz)
         type(interpolant) :: i

         i%siz = siz
         allocate( i%x(siz), i%y(siz) )
         i%x(:) = x(:) ; i%y(:) = y(:)
   end subroutine create_1d_interpolant

   !> evaluate interpolant at point xwant. at the moment only linear interpolation is supported
   !> but higher order methods could easily be implemented
   subroutine interpolate_1d( i, xwant_in, yrslt )
         type(interpolant) :: i !< interpolant
         real(r8) :: &
               xwant_in, &   !< point at which to evaluate interpolating function
               xwant, &
               yrslt, &      !< result (value of interpolating function at xwant)
               y1, y2, x1, x2
         integer :: idx

         ! clip xwant to table edgess
         xwant = max(xwant_in, i%x(1)    )
         xwant = min(xwant   , i%x(i%siz))

         ! perform bisection search for xwant in interpolant
         call bsearch_r8( i%x, xwant, idx )

         ! linear interpolation
         x1 = i%x(idx) ; x2 = i%x(idx+1)
         y1 = i%y(idx) ; y2 = i%y(idx+1)

         yrslt = y1 + ( y2 - y1 ) * ( xwant - x1 ) / ( x2 - x1 )
   end subroutine interpolate_1d


!   subroutine generic_interp_2d(xdim, ydim, xco, yco, n, f, myx, myy, ix1, ix2, iy1, iy2, res, dresdx, dresdy)
!
!         ! interpolate n 2d arrays of type real at the same location; return result and derivatives.
!         !
!         ! this routine is for multiple arrays because the interpolants only need to be calculated once per set of coordinates (e.g.
!         ! [rho,T]), otherwise there is a significant speed hit on the code
!         !
!         ! TODO: even better still may be to just have a routine to calculate the interpolants and then pass them to this routine
!         !
!         !     xdim, ydim           :  dimensions of data
!         !     xco, yco             :  x and y coordinates
!         !     f                    :  2d data of shape (n,xdim,ydim)
!         !     n                    :  number of 2d arrays to interpolate
!         !     myx, myy             :  x and y values for which result is desired
!         !     ix1, ix2, iy1, iy2   :  indices above/below desired point (this is not calculated automatically because different
!         !                             data types have different optimal methods for finding this, e.g. linear search vs bisection
!         !                             vs uniform spacing)
!         !     res, dresdx, dresdy  :  result and its derivatives wrt both x and for each interpolated array  (n)
!
!         integer, intent(in) :: xdim, ydim, ix1, ix2, iy1, iy2, n
!         real(r8), intent(in) :: xco(xdim), yco(ydim), f(n,xdim,ydim), myx, myy
!         real(r8) :: x1, x2, y1, y2, dx, dy, xa, xb, ya, yb
!         real(r8), dimension(n) :: f11, f21, f12, f22, fA, fB, fC, fD
!         real(r8), intent(out), dimension(n) :: res, dresdx, dresdy
!         integer :: i
!
!         x1 = xco(ix1); x2 = xco(ix2)
!         y1 = yco(iy1); y2 = yco(iy2)
!         dx = x2 - x1; dy = y2 - y1
!
!         xb = ( myx - x1 ) / dx; xa = 1._r8 - xb
!         yb = ( myy - y1 ) / dy; ya = 1._r8 - yb
!
!         f11(:) = f( :, ix1 , iy1 )
!         f21(:) = f( :, ix2 , iy1 )
!         f12(:) = f( :, ix1 , iy2 )
!         f22(:) = f( :, ix2 , iy2 )
!
!         fA(:)  = f11(:) * ya + f12(:) * yb
!         fB(:)  = f12(:) * xa + f22(:) * xb
!         fC(:)  = f21(:) * ya + f22(:) * yb
!         fD(:)  = f11(:) * xa + f21(:) * xb
!
!         res   (:)   = fA(:) * xa + fC(:) * xb
!         dresdx(:)   = ( fC(:) - fA(:) ) / dx
!         dresdy(:)   = ( fB(:) - fD(:) ) / dy
!
!   end subroutine generic_interp_2d
!

end module interp
