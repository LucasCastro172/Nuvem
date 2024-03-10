! Aluna: Amanda Guimarães Pereira
module integral_funcao
implicit none
Real:: R, p, L, a
Contains
!-----------------------------------
Real function funcao(x)
implicit none
Real, intent(in):: x
Real :: nu, d

nu =((p+L)**2+(a-x)**2)*(sqrt(R**2-x**2)+sqrt(p**2+R**2+ a**2-2*x*a))**2
d =(p**2+(a-x)**2)*(sqrt(R**2-x**2)+sqrt((p+L)**2+R**2+ a**2-2*x*a))**2
funcao=log(nu/d)

End Function
!-----------------------------------
Subroutine trapezios(F,a,b,e,S)
implicit none 
real, intent(in):: a,b,e
real,intent(out):: S ! saída
real,external:: F  ! quando a função está fora do modulo
integer:: i, N
real:: h,x
real:: soma,S1,S2,f_ab,Sn,S2n,S4n ! s2 interação atual, s1 interação anterior 
! S1
f_ab=f(a)+f(b)
N=4
h=(b-a)/N
x= a 
Soma=0.
Do i=1,N-1 
    x = x + h
    Soma = Soma + f(x)
End Do  
Sn = h*( 2*Soma + f_ab )/2.
! S2n
 f_ab=f(a)+f(b)
   N=8
   h=(b-a)/N
    x= a 
    Soma=0.
    Do i=1,N-1 
       x = x + h
       Soma = Soma + f(x)
    End Do   
    !S1 = 2*Soma
    S2n = h*( 2*Soma + f_ab )/2.    
    S1=(4*S2n-Sn)/3

Do 
    x= a + h/2. 
    Soma=Soma +f(x)
    Do i=1,N-1 
       x = x + h
       Soma = Soma + f(x)
    End Do   
    !h = h/2.
    S4n =h*( 2*Soma + f_ab )/4.
    !S2=(4*S4n-S2n)/3
    h = h/2.
    N = 2*N
    S2=(4*S4n-S2n)/3
    if(abs(s2-s1)>e) then
       S1=S2
       S2n=S4n
    Else
       S=S2
       exit
    End if           
End Do 
 End subroutine trapezios
!----------------------------------- 
End module

program integral
use integral_funcao
implicit none
Real, Allocatable,Dimension(:):: S
Real, Allocatable,Dimension(:):: x
Integer :: i,j
real:: eps
Integer:: N ! numeros de pontos

open(22,file='parametros.dat',status='old',action='read')
open(33,file='resultado_ex2.dat',status='replace',action='write') 

Read(22,*) R, L, p


Read(22,*) n 

Allocate(S(n),x(n))


eps=1.e-4

!x(1)=0
   Do j=1,n
     Read(22,*) a       
        x(j)= a
      Call trapezios(funcao,-R,R,eps,S(j))       
   End Do    

Do i=1,n
   write(*,*)x(i)
End do


Do i=1,n
   write(33,*)x(i),S(i)
End do

End program
