module integral 
implicit none 
Contains

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
write(*,*)'N:', N
 End subroutine trapezios
!!!!!!!!!!!!!!!!!!!!!!!
 real function parabola(x)
 implicit none
 real,intent(in):: x
  parabola=x**2
 End function
!!!!!!!!!!!!!!!!!!!!!!!
 real function seno(x)
 implicit none
 real,intent(in):: x 
 seno=sin(x)
 End function
 End module integral 
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!

program testa_trapezios
use integral
implicit none 
real:: a,b,s
real:: eps
real,parameter:: pi=3.14159265359

eps=1.e-5
!a = 0.
!b = 1.

a = 0.
b = pi/2
!N = 100

call trapezios(seno,a,b,eps,s)
write(*,*)'função y=x**2'
write(*,*) a,b
write(*,*)'Integral:',s
End program
