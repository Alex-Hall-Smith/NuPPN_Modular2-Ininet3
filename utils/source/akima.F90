! akima interpolation, based on routines from Andrew Ning:
! https://github.com/andrewning/akima
! 
! TODO replace linear search with binary search in "interp" subroutine!!
!
module akima
   use utils, only: i4, r8, bsearch_r8
   implicit none
   private

   type, public :: akima_interpolant
      integer(i4) :: n !< number of data points
      real(r8), allocatable, dimension(:) :: p0, p1, p2, p3, & !< polynomial coefficients
            xpt(:), ypt(:) !< x and y data points

   end type akima_interpolant

   public :: create_akima_interpolant, interp_akima

contains

   subroutine create_akima_interpolant(x, y, ai)
         implicit none
         real(r8) :: x(:), y(:)
         type(akima_interpolant) :: ai

         ai%n = size(x)
         allocate(ai%p0(ai%n-1), ai%p1(ai%n-1), ai%p2(ai%n-1), ai%p3(ai%n-1))
         allocate(ai%xpt(ai%n), ai%ypt(ai%n))

         ai%xpt = x; ai%ypt = y
         call setup(ai%n, ai%xpt, ai%ypt, ai%p0, ai%p1, ai%p2, ai%p3, 0.1_r8)
   end subroutine create_akima_interpolant

   ! evaluate an akima interpolant at a single x value (xx)
   subroutine interp_akima(ai, xx_in, res)
         implicit none
         type(akima_interpolant) :: ai
         real(r8) :: xx, xx_in, res, xpass(1), y(1), dydx(1), dydxpt(1), dydypt(1)

         ! clip input to table edges
         xx = max(xx_in, ai % xpt(1)   )
         xx = min(xx,    ai % xpt(ai%n))
         xpass(1) = xx
         call interp(ai%n, 1, xpass, ai%xpt, ai%p0, ai%p1, ai%p2, ai%p3, &
               y, dydx)

         res = y(1)

   end subroutine interp_akima

   subroutine abs_smooth(x, delta_x, y)
         ! absolute value function with small quadratic in the valley
         ! so that it is C1 continuous

         implicit none
         integer, parameter :: dp = kind(0.d0)

         ! in
         real(dp), intent(in) :: x, delta_x

         ! out
         real(dp), intent(out) :: y


         if (x >= delta_x) then
            y = x
         elseif (x <= -delta_x) then
            y = -x
         else
            y = x**2/(2.0_dp*delta_x) + delta_x/2.0_dp
         end if

   end subroutine abs_smooth



   subroutine setup(n, xpt, ypt, p0, p1, p2, p3, delta_x)
         ! setup the Akima spline function

         implicit none

         integer, parameter :: dp = kind(0.d0)
         real(dp), parameter :: eps = 1d-30

         ! in
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: xpt, ypt  ! given points
         real(dp), intent(in) :: delta_x
         !f2py real(dp), intent(in) :: delta_x = 0.1

         ! out
         real(dp), dimension(n-1), intent(out) :: p0, p1, p2, p3  ! spline coefficients

         ! local
         integer :: i
         real(dp), dimension(-1:n+1) :: m
         real(dp), dimension(n) :: t
         real(dp) :: m1, m2, m3, m4, w1, w2
         real(dp) :: t1, t2, dx

         ! compute segment slopes
         do i = 1, n-1
            m(i) = (ypt(i+1) - ypt(i)) / (xpt(i+1) - xpt(i))
         end do
!
         !! estimation for end points
         m(0) = 2.0_dp*m(1) - m(2)
         m(-1) = 2.0_dp*m(0) - m(1)
         m(n) = 2.0_dp*m(n-1) - m(n-2)
         m(n+1) = 2.0_dp*m(n) - m(n-1)

         ! slope at points
         do i = 1, n
            m1 = m(i-2)
            m2 = m(i-1)
            m3 = m(i)
            m4 = m(i+1)
            !         w1 = abs(m4 - m3)
            !         w2 = abs(m2 - m1)
            call abs_smooth(m4 - m3, delta_x, w1)
            call abs_smooth(m2 - m1, delta_x, w2)
            if ( w1 < eps .and. w2 < eps ) then
               t(i) = 0.5_dp*(m2 + m3)  ! special case to avoid divide by zero
            else
               t(i) = (w1*m2 + w2*m3) / (w1 + w2)
            end if
         end do

         ! polynomial cofficients
         do i = 1, n-1
            dx = xpt(i+1) - xpt(i)
            t1 = t(i)
            t2 = t(i+1)
            p0(i) = ypt(i)
            p1(i) = t1
            p2(i) = (3.0_dp*m(i) - 2.0_dp*t1 - t2)/dx
            p3(i) = (t1 + t2 - 2.0_dp*m(i))/dx**2
         end do

   end subroutine setup


   subroutine interp(npt, n, x, xpt, p0, p1, p2, p3, y, dydx)
         ! evaluate the Akima spline and its derivatives

         implicit none
         integer, parameter :: dp = kind(0.d0)

         ! in
         integer, intent(in) :: npt
         real(dp), dimension(npt), intent(in) :: xpt  ! given x points
         real(dp), dimension(npt-1), intent(in) :: p0, p1, p2, p3  ! spline coefficients
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x  ! x values to evalute at

         ! out
         real(dp), dimension(n), intent(out) :: y  ! interpolate y values
         real(dp), dimension(n), intent(out) :: dydx  ! derivative of y w.r.t. x

         ! local
         integer :: i, j, k
         real(dp) :: dx

         ! interpolate at each point
         do i = 1, n

            ! find location in array (use end segments if out of bounds)
            if (x(i) < xpt(1)) then
               j = 1

            else
               ! linear search for now
               do j = npt-1, 1, -1
                  if ( x(i) >= xpt(j)) then
                     exit
                  end if
               end do
            end if

            ! evaluate polynomial (and derivative)
            dx = (x(i) - xpt(j))
            y(i) = p0(j) + p1(j)*dx + p2(j)*dx**2 + p3(j)*dx**3
            dydx(i) = p1(j) + 2.0_dp*p2(j)*dx + 3.0_dp*p3(j)*dx**2

         end do


   end subroutine interp

end module akima
