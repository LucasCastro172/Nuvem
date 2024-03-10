!!!!!!!!!!!!!!!
Module Integrais
implicit none
Integer,parameter :: pr = kind(1.d0)
real(pr), parameter :: pi = 3.14159265359


contains
!-----------------------------------
Real(pr) function funcao(t)
implicit none
Real(pr), intent(in):: t
real (pr)           :: x
x = -(t*t)
funcao = 2/sqrt(pi)*exp(x)

End Function
!-----------------------------------
subroutine trapezios(FUN,a,b,S )
implicit none
real(pr),external       :: Fun
real(pr),intent(in)     :: a,b  !limites
real(pr),intent(out)    ::  S !resultado
real(pr)                :: soma,sex,h,x 
real(pr)                :: S1,S2
real(pr)                :: eps
integer                 :: i,N,k
        
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
        s2 = h * (sex + 2*soma )/2.d0
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
!-----------------------------------
End module
!!!!!!!!!!!!!!
program Funcao_erro
use integrais
implicit none
Real(pr)    :: x(60)       !eixo x
real(pr)    :: erro(60)    !resposta
Integer :: i
Real(pr)    :: xi

open(33,file='resultado_ex7.dat',status='replace',action='write') 


x(1) = -2.9
xi= -3
Do i=1,59
	X(i+1)=X(i) + 0.1
End Do

!write(*,*) x(60)

Do i=1,60  
	Call trapezios(Funcao,0.d0,x(i),erro(i))  
	write(*,*) erro(i)  
End Do  

Do i=1,60
	write(33,*) x(i), erro(i) 
End do


end program
