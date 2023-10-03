module cubinterp
   implicit none
   private

   ! ^_^ kinds a la IEEE754 via Reinecke
   integer, parameter, public :: &
         i4 = selected_int_kind  (9), &      !< four byte signed integer kind
         i8 = selected_int_kind  (15), &      !< eight byte signed integer kind
         r4 = selected_real_kind ( 6, 37), & !< four byte real kind (ieee754)
         r8 = selected_real_kind (15,307)    !< eight byte real kind (ieee754)

   real(r8), parameter, public :: &
         ZERO  = 0._r8, &
         ONE   = 1._r8, &
         SIXTH = ONE / 6._r8, &
         THIRD = ONE / 3._r8

   type, public :: interpolant_1d
      real(r8), allocatable :: &
            x(:), &   !< independent variable
            f(:), &   !< dependent variable
            fpp(:)    !< second derivative of independent variable wrt x
      integer(i4) :: &
            n         !< number of knots (x,f values)
   end type interpolant_1d

   type, public :: interpolant_2d
      real(r8), allocatable :: &
            x(:), &   !< independent variable
            y(:), &   !< independent variable
            f(:,:), & !< dependent variable
            fpp(:,:)  !< second derivative of independent variable in wrt x
      integer(i4) :: &
            nx, & !< number of x-knots
            ny    !< number of y-knots
   end type interpolant_2d

   public :: interp1d_cub, create_1d_cub_interpolant, interp2d_cub, create_2d_cub_interpolant

contains

   !> evaluate interpolating function at xwant
   subroutine evaluate_1d_cub_interpolant(x, f, fpp, xwant, res)
         use utils, only: bsearch_r8
         implicit none
         real(r8) :: x(:), f(:), fpp(:), xwant
         real(r8) :: res, xin, xl, xh, dx, A, B, C, D, Ap3, Bp3, dxp2
         integer(i4) :: il, ih

         ! schnips
         xin = clip(xwant, x)

         ! locate x in array
         call bsearch_r8( x, xin, il )
         ih = il + 1
         
         ! bounding x values
         xl = x(il); xh = x(ih)
         dx = xh - xl

         ! coefficients of cubic spline
         A = (xh - xin)/dx; B = 1. - A
         Ap3 = A**3; Bp3 = B**3; dxp2 = dx**2
         C = SIXTH * (Ap3 - A) * dxp2
         D = SIXTH * (Bp3 - B) * dxp2

         res = A * f(il) + B * f(ih) + C * fpp(il) + D * fpp(ih)
   end subroutine evaluate_1d_cub_interpolant

   !> perform 1d interpolation given 1d interpolant t and required coordinate xin
   !> to give result res
   subroutine interp1d_cub(t, xin, res)
         implicit none
         type(interpolant_1d) :: t
         real(r8) :: xin, res

         call evaluate_1d_cub_interpolant(t%x, t%f, t%fpp, xin, res)
   end subroutine interp1d_cub

   !> perform 2d interpolation given 2d interpolant t and required coordinates xin
   !> and yin to give result res
   subroutine interp2d_cub(t, xin, y, res)
         use utils, only: bsearch_r8
         implicit none
         type(interpolant_2d) :: t
         type(interpolant_1d) :: tcol
         real(r8) :: xin, y, yin, res, ynew(5), fnew(5)
         integer(i4) :: j, iy, i
         
         ! find and y
         yin = clip(y, t%y)
         call bsearch_r8(t%y, yin, iy)

         ! create new column of data with lenth 5 by cubic interpolation in x and construct an
         ! interpolant in y
         iy = max(3, min(iy,t%ny-2)); i = 0
         do j = iy-2, iy+2
            i = i + 1
            call evaluate_1d_interpolant(t%x, t%f(:,j), t%fpp(:,j), xin, fnew(i))
         end do
         ynew(:) = t%y(iy-2:iy+2)
         call create_1d_interpolant(tcol, ynew, fnew, ZERO, ZERO )

         ! evaluate new interpolant at required y
         call interp1d_cub(tcol, yin, res)
   end subroutine interp2d_cub



   !> solve the linear system for the cubic spline interpolation problem
   subroutine calculate_1d_cub_interpolant(x, f, fpp, fpp1, fppn)
         implicit none
         integer(i4) :: i, s
         real(r8) :: x(:), f(:), fpp(:), &
               fpp1, & !< boundary contition 1: second derivative of f wrt x at 1: 0 = natural
               fppn    !< boundary contition 2: second derivative of f wrt x at n: 0 = natural
         real(r8), allocatable :: a(:), b(:), g(:), rhs(:)

         s = size(x)
         ! allocate storage for rhs and banded (tridiagonal) matrix
         allocate(rhs(2:s-1), a(2:s-1), b(2:s-1), g(2:s-1))

         ! calculate rhs and xtri
         do i = 2, s-1
            a(i) = SIXTH * (x(i) - x(i-1))
            b(i) = THIRD * (x(i+1) - x(i-1))
            g(i) = SIXTH * (x(i+1) - x(i))

            rhs(i) = (f(i+1) - f(i)) / (x(i+1) - x(i)) - (f(i) - f(i-1)) / (x(i) - x(i-1))
         end do

         ! boundary conditions
         fpp(1) = fpp1; fpp(s) = fppn

         ! solve the system for ypp
         call solve_tridiag(s-2, a, b, g, rhs, fpp(2:s-1))
   end subroutine calculate_1d_cub_interpolant


   subroutine create_1d_cub_interpolant(t, x, f, fpp1, fppn)
         implicit none
         type(interpolant_1d) :: t
         real(r8), intent(in) :: &
               x(:), & !< independent variable
               f(:), & !< dependent variable
               fpp1, & !< boundary contition 1: second derivative of f wrt x at 1: 0 = natural
               fppn    !< boundary contition 2: second derivative of f wrt x at n: 0 = natural
         integer(i4) :: s

         ! store x and f
         s = size(x); t%n = s
         allocate(t%x(s), t%f(s), t%fpp(s))
         t%x = x; t%f = f

         call calculate_1d_cub_interpolant(t%x, t%f, t%fpp, fpp1, fppn)
   end subroutine create_1d_cub_interpolant


   subroutine create_2d_cub_interpolant(t, x, y, f)
         implicit none
         type(interpolant_2d) :: t
         real(r8), intent(in) :: x(:), y(:), f(:,:)
         integer(i4) :: sx, sy, i

         ! store x, y and f
         sx = size(x); sy = size(y); t%nx = sx; t%ny = sy
         allocate(t%x(sx), t%y(sy), t%f(sx,sy), t%fpp(sx,sy))
         t%x = x; t%y = y; t%f = f

         do i = 1, sy
            call calculate_1d_cub_interpolant(t%x, t%f(:,i), t%fpp(:,i), ZERO, ZERO)
         end do
   end subroutine create_2d_cub_interpolant


   !> clip a to the bounds of array b
   real(r8) function clip(a,b)
         real(r8) :: a, b(:)
         integer(i4) :: n

         n = size(b)
         if (b(1) < b(n)) then
            clip = min(max(a, b(1)), b(n))
         else
            clip = min(max(a, b(n)), b(1))
         end if
   end function clip
   

   !> solve tridiagonal system. a, b, c are the diagonals, x will be the solution
   subroutine solve_tridiag(n, a, b, c, rhs, x)
         integer(i4) :: n, i
         real(r8) :: a(n), b(n), c(n), rhs(n), x(n), inv, dx

         ! forward elimination
         c(1) = c(1) / b(1)
         rhs(1) = rhs(1) / b(1)
         do i = 2, n
            inv  = ONE / (b(i) - a(i) * c(i-1))
            c(i) = c(i) * inv
            rhs(i) = (rhs(i) - a(i) * rhs(i-1) ) * inv
         end do

         ! back substitution
         x(n) = rhs(n)
         do i = n-1, 1, -1
            x(i) = rhs(i) - c(i) * x(i+1)
         end do
   end subroutine solve_tridiag


end module cubinterp
