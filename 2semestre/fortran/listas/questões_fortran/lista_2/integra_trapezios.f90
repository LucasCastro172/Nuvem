! integração numérica pelo método dos trapézios 
module integral 
implicit none 

Contains

Subroutine trapezios(F,a,b,N,s)
implicit none 
real, intent(in):: a,b
integer, intent(in):: N
real,intent(out):: S ! saída
real,external:: F  ! quando a função está fora do modulo
integer:: i
real:: h,x

h=(b-a)/N

!x= a + h
!S=0.
!  Do i=1,N-1
!     S = S + f(x)
!     x = x + h
!  End Do   

 x= a 
 S=0.
  Do i=1,N-1 
     x = x + h
     S = S + f(x)
  End Do   
 S = 2*S
 S = h*( f(a) + S + f(b) )/2
 
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
integer:: N
real,parameter:: pi=3.14159265359

a = 0.
b = 2*pi
N = 100

call trapezios(seno,a,b,N,s)
write(*,*)'função seno(x)'
write(*,'("Intervalo:",2(f15.5))') a,b
write(*,*)'N:', N
write(*,*)'Integral:',s


End program
