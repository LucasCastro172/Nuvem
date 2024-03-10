!!!!!!!!!!!!!!
Module sistemas
implicit none

integer, parameter :: prc = kind(0.d0)

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


subroutine elimina(n,a,b)
!essa versão n está verificando a qualidade da matriz
!pressupoe q é uma matriz singular
implicit none
integer,  intent(in)     :: n
real(prc), intent(inout) :: A(n,n), b(n)
real(prc)                :: m
integer                  :: i,j,k

do k = 1, N-1
	do i = k+1, n
		m = a(i,k) / a(k,k)
		!a(i,:) = a(i,:) - m * A(k,:)
		do j = k+1,n
			a(i,j) = a(i,j) - m * A(k,j)
		end do 
		b(i) = b(i) - m*b(k)
	end do
end do	



end subroutine  elimina


end module


!!!!!!!!
!!!!!!!

Program teste_sistema
use sistemas
implicit none
integer, parameter :: n=3
real(8)            :: a(n,n), b(n), x(n), x0(n)
integer            :: i,j

A(1,:) = [2.d0, 1.d0, -1.d0]
A(2,:) = [-1.d0, 4.d0, 1.d0]
A(3,:) = [1.d0,1.d0,1.d0]
b = [3.d0,5.d0,8.d0]

call elimina(n,a,b)

do i = 1,n
	write(*,'(4())')
end do

end program
