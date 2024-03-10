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


subroutine elimina_pivoteamento(n,a,b,erro)
!essa versão n está verificando a qualidade da matriz
!pressupoe q é uma matriz singular
implicit none
integer,  intent(in)     :: n
real(prc), intent(inout) :: A(n,n), b(n)
real(prc)                :: m
integer,  intent(out)    :: !  Se erro igual 
real(prc)                :: eps =  1.d-8
real(prc)                :: aux, pivo
integer                  :: i,j,k, lp

do k = 1, N-1
!pivoteamento
	erro = 0
	pivo = a(k,k)
	lp = k
	do i = k+1,n
		if (abs(a(i,k)) > abs(pivo) ) then
			pivo = a(i,k)
			lp = i
		end if
	end do
	if ( abs(pivo) > eps)  then
		if (lp /= k) then
			do  j  = k,n
				aux = a(k,K)
				a(k,k) = a(lp,j)
				a(lp,j) =  aux
			end do
			aux = b(k)
			b(k) = b(lp)
			b(lp) = aux
		end if
!!!!!!!!!!!!!!!!!!!Eliminação
		do i = k+1, n
			m = a(i,k) / a(k,k)
		!a(i,:) = a(i,:) - m * A(k,:)
			do j = k+1,n
				a(i,j) = a(i,j) - m * A(k,j)
			end do 
			b(i) = b(i) - m*b(k)
			end do
		end do		
	else 
	erro = 1
	end if
end do

if (abs(a(n,n)) < eps) erro =  1

end subroutine  elimina_pivoteamento


subroutine fatora_U(n, a)
!pressupoe uma matriz u do tipo triangular superior
implicit none
integer, intent(in)   :: n
real(prc), intent(inout) :: a(n,n)
integer               :: i,j,k
real(prc)             :: m


do k = 1, n-1
	do i = k+1, n
		m = a(i,k) / a(k,k)
		do j = k+1,n
			a(i,j) = a(i,j) - m * A(k,j)
		end do 
		a(i,k) = m
	end do
end do		

end subroutine fatora_u

subroutine subst_direta(n,l,b,y)
implicit none
integer, intent(in)    :: n
real(prc), intent(in)  :: L(n,n), b(n)
real(prc), intent(out) :: y(n)
integer                :: i,j
real(prc)              :: soma

y(1) = b(1)

do i=2, n
	soma = 0.d0
	do j = 1, i-1
		soma = soma + l(i,j)*y(j)
	end do
	y(i) = 
end do


end subroutine subst_direta


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
A(3,:) = [1.d0,1.d0,0.d0]
b = [3.d0,5.d0,8.d0]

call elimina_pioteamento(n,a,b,erro)
write(*,*) 

do i = 1,n
	write(*,'(4())')
end do

end program
