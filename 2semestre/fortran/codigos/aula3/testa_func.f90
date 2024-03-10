Real function exp_cos (x)
implicit none
real, intent(In) :: x

exp_cos= exp(-x**2) * cos(10.*x)


end function exp_cos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program testafunc
implicit none
real, external :: exp_cos
real :: a

write(*,*) 'Valor teste'
read(*,*) a

write(*,*) 'Funcao:', exp_cos(a)


end program
