module constants
integer, parameter :: N = 128
integer, parameter :: M = 48

integer, parameter :: cx = N / 4
integer, parameter :: cy = M / 2
integer, parameter :: r = M / 9

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

integer, dimension(9) :: noslip = (/0, 2, 1, 6, 8, 7, 3, 5, 4/)

end module

subroutine update_rho(rho, f)
   use constants
   implicit none
   real(kind=8), dimension(N, M) :: rho
   real(kind=8), dimension(N, M, 9) :: f
   integer :: i, j, q
   real(kind=8) :: s

   !$omp parallel do collapse(2) private(i,j,q,s)
   do j=1,M
      do i=1,N
         s = 0.0d0
         do q=1,9
            s = s + f(i, j, q)
         end do

         rho(i, j) = s

      end do
   end do

end subroutine

subroutine update_u(u, rho, f)
   use constants
   implicit none
   real(kind=8), dimension(N, M, 2) :: u
   real(kind=8), dimension(N, M) :: rho
   real(kind=8), dimension(N, M, 9) :: f
   integer :: i, j, q
   real(kind=8) :: s1, s2

   !$omp parallel do collapse(2) private(i,j,q,s1,s2)
   do j=1,M
      do i=1,N
         s1 = 0.0d0
         s2 = 0.0d0
         do q=1,9
            s1 = s1 + f(i, j, q) * c(q, 1)
            s2 = s2 + f(i, j, q) * c(q, 2)
         end do

         u(i, j, 1) = s1 / rho(i, j)
         u(i, j, 2) = s2 / rho(i, j)

      end do
   end do

end subroutine

subroutine update_feq(feq, rho, u, cs)
   use constants
   implicit none
   real(kind=8), dimension(N, M, 9) :: feq
   real(kind=8), dimension(N, M) :: rho
   real(kind=8), dimension(N, M, 2) :: u
   real(kind=8) :: cs, oocs2, cu, u2
   integer :: i, j, q

   oocs2 = 1.0d0 / cs**2

   !$omp parallel do collapse(3) private(i,j,q,cu,u2)
   do q=1,9
      do j=1,M
         do i=1,N
            cu = c(q, 1) * u(i, j, 1) + c(q, 2) * u(i, j, 2)
            u2 = u(i, j, 1)**2 + u(i, j, 2)**2

            feq(i, j, q) = w(q) * rho(i, j) * (1.0d0 + oocs2*cu + 0.5d0*oocs2**2 * cu**2 - 0.5d0*oocs2 * u2)
         end do
      end do
   end do

end subroutine

subroutine collide(ft, f, feq, obstacle, omega)
   use constants
   implicit none
   real(kind=8), dimension(N, M, 9) :: ft, f, feq
   integer, dimension(N, M) :: obstacle
   real(kind=8) :: omega
   integer :: i, j, q

   !$omp parallel do collapse(3) private(i,j,q)
   do q=1,9
      do j=1,M
         do i=1,N
            if(obstacle(i,j) == 0) then
               ft(i,j,q) = f(i,j,q) - omega * (f(i,j,q) - feq(i,j,q))
            else
               ft(i,j,q) = f(i,j,noslip(q)+1)
            end if
         end do
      end do
   end do

end subroutine

subroutine stream(f, ft)
   use constants
   implicit none
   real(kind=8), dimension(0:N-1, 0:M-1, 9) :: f, ft
   integer :: i, j, q

   !$omp parallel do collapse(3) private(i,j,q)
   do q=1,9
      do j=0,M-1
         do i=0,N-1
            f(mod(N + mod(i + c(q,1), N), N), mod(M + mod(j + c(q,2), M), M), q) = ft(i, j, q)
         end do
      end do
   end do

end subroutine

function inlet_vel(k, i, ulb) result(v)
   use constants
   implicit none
   integer, intent(in) :: k, i
   real(kind=8), intent(in) :: ulb
   real(kind=8) :: x, v

   x = real(i,8) / (M-1)
   v = (1.0d0 - k) * ulb * (1.0d0 + 0.0001d0 * sin(2.0d0 * 3.142d0 * x))

end function

subroutine update(f, ft, feq, rho, u, obstacle, cs, omega, ulb)
   use constants
   implicit none
   real(kind=8), dimension(N, M, 9) :: f, ft, feq
   real(kind=8), dimension(N, M) :: rho
   integer, dimension(N, M) :: obstacle
   real(kind=8), dimension(N, M, 2) :: u
   real(kind=8) :: cs, omega, ulb, s1, s2
   integer :: j, q
   real(kind=8) :: inlet_vel

   ! outflow
   !$omp parallel do collapse(2) private(j,q)
   do q=4,6
      do j=1,M
         f(N, j, q) = f(N-1, j, q)
      end do
   end do

   call update_rho(rho, f)
   call update_u(u, rho, f)

   ! inflow
   !$omp parallel do private(j,q,s1,s2)
   do j=1,M
      s1 = 0.0d0
      s2 = 0.0d0
      do q=1,6
         if(q <= 3) then
            s2 = s2 + f(1, j, q)
         else
            s1 = s1 + f(1, j, q)
         end if
      end do

      u(1, j, 1) = inlet_vel(0, j, ulb)
      u(1, j, 2) = inlet_vel(1, j, ulb)

      rho(1, j) = 1.0d0 / (1.0d0 - u(1, j, 1)) * (s2 + 2.0d0 * s1)
   end do
   call update_feq(feq, rho, u, cs)
   !$omp parallel do collapse(2) private(j,q)
   do q=7,9
      do j=1,M
         f(1, j, q) = feq(1, j, q)
      end do
   end do

   call collide(ft, f, feq, obstacle, omega)
   call stream(f, ft)

end subroutine

subroutine print_u(u, t)
   use constants
   implicit none
   real(kind=8), dimension(N, M, 2) :: u
   integer :: t, i, j

   write (*, '(i0,x)', advance='no') t
   do i=1,N
      do j=1,M
         write (*, '(f0.10,x)', advance='no') sqrt(u(i,j,1)**2 + u(i,j,2)**2)
      end do
   end do
   print *

end subroutine

program lb
   use constants
   implicit none
   real(kind=8), parameter :: cs = 1.0d0 / sqrt(3.0d0)
   real(kind=8), parameter :: Re = 160.0d0
   real(kind=8), parameter :: ulb = 0.04d0
   real(kind=8), parameter :: nu = ulb * r / Re
   real(kind=8), parameter :: omega = 1.0d0 / (3.0d0 * nu + 0.5d0)
   integer, parameter :: T = 20000

   real(kind=8), allocatable :: f(:,:,:), ft(:,:,:), feq(:,:,:)
   real(kind=8), allocatable :: rho(:,:)
   integer, allocatable :: obstacle(:,:)
   real(kind=8), allocatable :: u(:,:,:)

   integer :: i, j, q, time
   real(kind=8) :: inlet_vel

   allocate(f(N,M,9))
   allocate(ft(N,M,9))
   allocate(feq(N,M,9))
   allocate(rho(N,M))
   allocate(obstacle(N,M))
   allocate(u(N,M,2))

   do j=1,M
      do i=1,N
         rho(i,j) = 1.0d0
         u(i,j,1) = inlet_vel(0, j, ulb)
         u(i,j,2) = inlet_vel(1, j, ulb)
         if((i-cx)**2 + (j-cy)**2 < r**2) then
            obstacle(i,j) = 1
         else
            obstacle(i,j) = 0
         end if
      end do
   end do

   call update_feq(feq, rho, u, cs)

   do q=1,9
      do j=1,M
         do i=1,N
            f(i,j,q) = feq(i,j,q)
         end do
      end do
   end do

   do time=0,T-1
      if(mod(time, 100) == 0) then
         call print_u(u, time)
      end if
      call update(f, ft, feq, rho, u, obstacle, cs, omega, ulb)
   end do

end program
