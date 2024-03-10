!!!!!!!!!!
program equacao
implicit none
real :: a, b, c
real :: r1, r2
integer :: nr
real :: degree, r
real, parameter :: pi = 3.14159265359




degree = r * (180/pi)

write(*,*) 'calculo das raizes'
write(*,*) 'Valores A, B e C'
read(*,*) a, b, c





call eq2grau(a, b, c, nr, r1, r2)

if (nr == 2) then
	write(*,*) 'Duas raizes:'
	write(*,*) 'r1:', r1
	write(*,*) 'r2:', r2
else if (nr == 1) then
	write(*,*) 'Uma raiz:'
	write(*,*) 'r1:', r1
else if (nr  == 0) then
	write(*,*) 'não há raizes reais'
end if


CONTAINS

subroutine eq2grau( a, b, c, nr, r1, r2)
!calculo das raizes da equacao ax^2 + bx + c =0
implicit none
real, intent(in) :: a, b, c
integer, intent(out) :: nr
real, intent(out) ::r1, r2
real :: delta

Delta = b**2 - 4.*a*c

if (Delta > 0.) then
	nr = 2
	r1 =  (-b + sqrt(delta))/(2.*a)
	r2 =  (-b - sqrt(delta))/(2.*a)
else if (delta < 0.) then
	nr = 0
else
	nr = 1
	r1 = -b/(2.*a)
	r2 = r1
end if

end subroutine eq2grau
	
	
	
	
	
	



end program	
