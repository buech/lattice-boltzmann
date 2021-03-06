module constants
real, parameter :: cx = N / 4.0
real, parameter :: cy = M / 2.0
real, parameter :: r = M / 9.0

real(kind=8), dimension(9) :: w = (/&
   4.0d0/9.0d0 ,&
   1.0d0/9.0d0 ,&
   1.0d0/9.0d0 ,&
   1.0d0/9.0d0 ,&
   1.0d0/36.0d0,&
   1.0d0/36.0d0,&
   1.0d0/9.0d0 ,&
   1.0d0/36.0d0,&
   1.0d0/36.0d0 &
   /)

integer, dimension(9,2) :: c = transpose(reshape((/&
   0,  0,&
   0, -1,&
   0,  1,&
  -1,  0,&
  -1, -1,&
  -1,  1,&
   1,  0,&
   1, -1,&
   1,  1 &
   /), (/2, 9/)))

integer, dimension(9) :: noslip = (/1, 3, 2, 7, 9, 8, 4, 6, 5/)

end module

function pmod(a, b) result(res)
   implicit none
   integer, intent(in) :: a, b
   integer :: res

   res = mod(b + mod(a,b), b)

end function

function inlet_vel(i) result(v)
   implicit none
   integer :: i
   real(kind=8) :: x, v

   x = real(i,8) / (M-1)
   v = ULB * (1 + 0.0001d0 * sin(2 * 3.142d0 * x))

end function

function rho(i, j, f) result(s)
   implicit none
   real(kind=8), dimension(9, N, M) :: f
   integer :: i, j, q
   real(kind=8) :: s

   s = 0.0d0
   do q=1,9
      s = s + f(q, i, j)
   end do

end function

function u(k, i, j, f) result(s)
   use constants, only : c
   implicit none
   real(kind=8), dimension(9, N, M) :: f
   integer :: i, j, q, k
   real(kind=8) :: s, rho

   s = 0.0d0
   do q=1,9
      s = s + f(q, i, j) * c(q, k)
   end do

   s = s / rho(i, j, f)

end function

subroutine boundary(f)
   use constants, only : c, w
   implicit none
   real(kind=8), dimension(9, N, M) :: f
   real(kind=8) :: ux, s1, s2, cu, rho_j, inlet_vel
   integer :: j, q

   ! outflow
   do q=4,6
      do j=1,M
         f(q, N, j) = f(q, N-1, j)
      end do
   end do

   ! inflow
   do j=1,M
      ux = inlet_vel(j)
      s2 = f(1, 1, j) + f(2, 1, j) + f(3, 1, j)
      s1 = f(4, 1, j) + f(5, 1, j) + f(6, 1, j)

      rho_j = 1.0d0 / (1 - ux) * (s2 + 2*s1)

      do q=7,9
         cu = c(q,1) * ux
         f(q, 1, j) = w(q) * rho_j * (1 + 3*cu + 0.5*9*cu**2 - 0.5*3*ux**2)
      end do
   end do

end subroutine

subroutine collstream(fnew, fold, obstacle, omega)
   use constants, only : c, w, noslip
   implicit none
   real(kind=8), dimension(9, N, M) :: fnew, fold
   integer, dimension(N, M) :: obstacle
   real(kind=8) :: omega, ux, uy, rho_ij, u2, cu
   integer :: i, j, q, pmod

   !$omp parallel do collapse(2) private(i,j,q,rho_ij,ux,uy,u2,cu)
   do j=1,M
      do i=1,N
         rho_ij = 0.0d0
         ux = 0.0d0
         uy = 0.0d0
         do q=1,9
            rho_ij = rho_ij + fold(q,i,j)
            ux = ux + fold(q,i,j) * c(q,1)
            uy = uy + fold(q,i,j) * c(q,2)
         end do
         ux = ux / rho_ij
         uy = uy / rho_ij
         u2 = ux**2 + uy**2
         do q=1,9
            cu = ux*c(q,1) + uy*c(q,2)
            if(obstacle(i,j) == 0) then
               fnew(q, pmod(i + c(q,1) - 1, N) + 1, pmod(j + c(q,2) - 1, M) + 1) =&
                  (1-omega) * fold(q,i,j)&
                  + omega * w(q) * rho_ij * (1 + 3*cu + 0.5*9*cu**2 - 0.5*3*u2)
            else
               fnew(q, pmod(i + c(q,1) - 1, N) + 1, pmod(j + c(q,2) - 1, M) + 1) = fold(noslip(q), i, j)
            end if
         end do
      end do
   end do

end subroutine

subroutine update(fnew, fold, obstacle, omega)
   implicit none
   real(kind=8), dimension(9, N, M) :: fnew, fold
   integer, dimension(N, M) :: obstacle
   real(kind=8) :: omega

   call collstream(fnew, fold, obstacle, omega)
   call boundary(fnew)

end subroutine

subroutine print_u(f, t)
   implicit none
   real(kind=8), dimension(9, N, M) :: f
   integer :: t, i, j
   real(kind=8) :: u

   write (*, '(i0,1x)', advance='no') t
   do i=1,N
      do j=1,M
         write (*, '(f0.10,1x)', advance='no') sqrt(u(1,i,j, f)**2 + u(2,i,j, f)**2)
      end do
   end do
   print *

end subroutine

program lb
   use constants, only : cx, cy, r, c, w
   implicit none
   real(kind=8), parameter :: nu = ULB * r / RE
   real(kind=8), parameter :: omega = 1.0d0 / (3.0d0 * nu + 0.5d0)

   real(kind=8), allocatable :: fnew(:,:,:), fold(:,:,:)
   integer, allocatable :: obstacle(:,:)

   integer :: i, j, q, time
   real(kind=8) :: cu, ux, inlet_vel

   allocate(fnew(9,N,M))
   allocate(fold(9,N,M))
   allocate(obstacle(N,M))

   do j=1,M
      do i=1,N
         if((i-cx)**2 + (j-cy)**2 < r**2) then
            obstacle(i,j) = 1
         else
            obstacle(i,j) = 0
         end if
      end do
   end do

   do j=1,M
      do i=1,N
         ux = inlet_vel(j)
         do q=1,9
            cu = ux * c(q,1)
            fold(q, i, j) = w(q) * (1 + 3*cu + 0.5*9*cu**2 - 0.5*3*ux**2)
         end do
      end do
   end do

   do time=0,T-1
      if(mod(time, 2) == 1) then
         call update(fold, fnew, obstacle, omega)
      else
         if(mod(time, 100) == 0) then
            call print_u(fold, time)
         end if
         call update(fnew, fold, obstacle, omega)
      end if
   end do

end program
