program exponencial
implicit none
integer :: i
real :: x, t 
real :: soma
real :: eps
integer :: N


write(*,*) 'x?'
read(*,*) x


eps = 1.e-5
soma = 1.
n = 0 
do 
	n = n +1
	t = x**N / fatorial(N)
	soma = soma + t
	write(*,*) n, t, soma
	if (abs(t) < eps) exit
end do

write(*,*) '  soma:', soma
write(*,*) 'exp(x):', exp(x)
	

Contains

integer function fatorial (N)
implicit none
integer, intent(in) :: N
integer :: i

fatorial = 1
do i = 2, N
	fatorial = fatorial * i
end do

end function fatorial



end program
