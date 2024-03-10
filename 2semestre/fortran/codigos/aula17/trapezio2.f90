!!!!!!!!!!!!!!!
Module Integrais
implicit none
Integer,parameter :: pr = kind(1.d0)
real(pr), parameter :: pi = 3.14159265359


contains

subroutine trapezios( FUN, a, b, s)
implicit none
real(pr), external :: fun
real(pr), intent(in):: a, b
real(pr), intent(out) :: s
real(pr) ::  soma, sex, h,x, s1, eps
integer :: i, n

eps = 1.d-6
n = 4
h = (b-a) / N

soma  = 0.
x = a + h
do i  = 1, N-1
	soma = soma + Fun(x)
	x = x + h
end do

s1 = h * ( fun(a) +  2*soma + fun(b) ) /  2.d0


do 
	x = a - h/2.d0
	do i = 1, N
		x = x  + h
		soma = soma + Fun(x)
	end do
	h = h / 2.d0
	s = h * (sex + 2*soma) / 2.d0
	
	if (abs(s - s1) > eps ) then
		n = 2*n
		s1 = s
	else
		exit
	end if
	

end do



end subroutine trapezios

!!!!!!!!!!!!
real(pr) function x_quad(x)
implicit none
real(pr), intent(in) :: x

x_quad = x**2

end function x_quad

real(pr) function seno(x)
implicit none
real(pr), intent(in) :: x

seno = sin(x)

end function seno
!!!!!!!!!!!!!
!real(pr) function x_inv(x)
!implicit none
!real(pr), intent(in) :: x

!x_inv = 

!end function x_inv


end module integrais


!!!!!!
Program Teste_trapezios
use integrais
implicit none
real(pr)  :: a,b
integer :: N
real(pr) :: area, x

!a = 0.
!b = 1.
!n =  100000



CALL TRAPEZIOS(X_QUAD,A, B, N, AREA)

WRITE(*,*) 'RESULTADO PARA X**2:'
WRITE(*,*) 'VALOR VERDADEIRO:',  1./3.
WRITE(*,*) 'VALOR ESTIMADO:', AREA

A = 0.
B = PI
N =  100000

CALL TRAPEZIOS(SENO,A, B, N, AREA)

WRITE(*,*) 'RESULTADO PARA SENO:'
WRITE(*,*) 'VALOR VERDADEIRO:',  2.
WRITE(*,*) 'VALOR ESTIMADO:', AREA







end program




