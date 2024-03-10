!!!!!!!!!!!!!!!
Module Integrais
implicit none
Integer,parameter :: pr = kind(1.d0)
real(pr), parameter :: pi = 3.14159265359


contains

subroutine trapezios(FUN,a,b,S )
        implicit none
        real(pr),external       :: FUN
        real(pr),intent(in)     :: a,b
        real(pr),intent(out)    ::  S
        real(pr)                :: soma,sex,h,x
        real(pr)                :: S1,S2
        real(pr)                :: eps
        integer             :: i,N,k
        
        eps = 1.d-6
        N = 4
        h = (b-a)/N
        sex = FUN(a) + FUN(b)
        soma = 0.
        x = a+h
        do i = 1, N-1
            soma = soma + FUN(x)
            x = x + h
        end do
        S1 = h * (sex + 2*soma)/2.d0
        
        k = 0
        do 
            x = a - h/2.d0
            do i = 1, N
                x = x + h
                soma = soma + FUN(x)
            end do
            h = h / 2.d0
            S2 = h * (sex + 2*soma )/2.d0
            S = (1.d0/3.d0) * (4*S2 - S1)
            if (abs(S-S1) > eps) then
                N = 2 * N
                S1 = S
            else
                exit
            end if
            k = k + 1
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


end module integrais


!!!!!!
Program Teste_trapezios
use integrais
implicit none
real(pr)  :: a,b
integer :: N
real(pr) :: area, x

a = 0.
b = 1.
!n =  100000



call trapezios(x_quad,a, b, area)

write(*,*) 'Resultado para x**2:'
write(*,*) 'Valor verdadeiro:',  1./3.
write(*,*) 'Valor estimado:', area

a = 0.
b = pi
n =  100000

call trapezios(seno,a, b, area)

write(*,*) 'Resultado para seno:'
write(*,*) 'Valor verdadeiro:',  2.
write(*,*) 'Valor estimado:', area







end program




