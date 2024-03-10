!!!!!!!!!!!!!!
Module sistemas
implicit none

contains

subroutine resolve_U(n, u, b, x)
!pressupoe uma matriz u do tipo triangular superior
implicit none
integer, intent(in) :: n
real(8), intent(in) :: u(n,n), b(n)
real(8), intent(out):: x(n)
integer             :: i,j
real(8)             :: soma

x(n) = b(n) / u(n,n)
do i = n-1, 1, -1
	soma = 0.d0
	do j = i+1, n
		soma = soma  + u(i,j) * x(j)
	end do
	x(i) = (b(i) - soma) / u(i,i)
end do


end subroutine resolve_u
end module


!!!!!!!!
!!!!!!!

Program teste_sistema
use sistemas
implicit none
integer, parameter :: n=3
real(8)            :: a(n,n), b(n), x(n), x0(n)
integer            :: i,j

!a = 0.d0
!A(1,1) = 2.d0
!A(1,2) = 1.d0
!A(1,3) = -1.d0
!A(2,2) = 9.d0
!a(2,3) = 1.d0
!a(3,3) = 13.d0
!b = [3.d0,]

!a(1,:) = [1.d0, 2.d0, 4.d0]
!a(2,:) = [0.d0, 3.d0, 5.d0]
!a(3,:) = [0.d0, 0.d0, 6.d0]

a = 0.d0

do j = 1,n 
	call random_number(A(1:j,j))
end do 
x0 = 1.d0

do i = 1, n 
	write(*,)


b = matmul(a,x0)

call resolve_u(n,a,b,x)


write(*, '(*(f15.8))') x 








end program 

