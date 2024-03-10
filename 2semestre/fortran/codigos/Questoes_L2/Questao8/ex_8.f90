module integral_funcao
implicit none
Real:: A, B, C
Contains
!-----------------------------------
Real function funcao(l)
Real, intent(in):: l
integer::j

funcao=A*exp(-B*l)*(B*cos(C*l)+C*sin(C*l))
!funcao=A(j)*exp(-B(j)*l)*(B(j)*cos(C(j)*l)+C(j)*sin(C(j)*l))

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
!!!!!!!!!!!!!!
program integral
use integral_funcao
implicit none
Real, Allocatable,Dimension(:,:):: S
Real, Allocatable,Dimension(:):: x
Real:: xi
Integer :: i,j
real:: eps
Allocate(S(50,7), x(50))

open(22,file='parametros_ABC.dat',status='old',action='read')
open(33,file='resultado_ex8.dat',status='replace',action='write') 

eps=1.e-5
x(1)=0.1
xi=0.
Do i=1,49
  X(i+1)=X(i) + 0.1
End Do
 
 Do j=1,7
    !Read(22,*) A(j), B(j), C(j)
    Read(22,*) A, B, C
      ! write(*,*) A,B,C
    Do i=1,50  
      Call trapezios(Funcao,xi,x(i),eps,S(i,j))    
    End Do  
 End Do    
 

  Do i=1,50
     write(33,*)x(i),(S(i,j),j=1,7)
  End do
End program
