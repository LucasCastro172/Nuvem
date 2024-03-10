Real function seno(x)
implicit none 

Real, intent(inout):: x
Integer:: N,i
Real:: soma
Real:: t
Real, parameter:: eps= 1.e-6
Real, parameter:: pi=3.14159265359 


If(2*pi<x) then
Do 
  x=x-2*pi
  if(x<0.) exit
  x=x-2*pi
 ! End if
  write(*,*)'x',x
  seno=x
  soma=x
  t=x
  N=1
  
  Do
     N=N+2
     t=-t*x**2/(N*(N-1))
     soma=soma+t
     If(abs(soma-seno)<eps) then
        seno=soma
        exit
     Else
        seno=soma
     End If
             
  End Do
End do
!!!!!!!!!!!!!!!!!!!11
!Else if(x<0) then
!Do 
!  x=x+2*pi
!  if(x>0.) exit
!  x=x+2*pi
! ! End if
!  write(*,*)'x',x
 ! seno=abs(x)
 ! soma=abs(x)
 ! t=abs(x)
 ! N=1
  
 ! Do
  !   N=N+2
   !  t=-t*x**2/(N*(N-1))
   !  soma=soma+t
   !  If(abs(soma-seno)<eps) then
   !     seno=-soma
   !     exit
   !  Else
   !     seno=-soma
   !  End If
             
  !End Do
!End do
!!!!!!!!!!!!!!!!!!!!!! 
Else
    seno=x
    soma=x
    t=x
    N=1  
    Do
     N=N+2
     t=-t*x**2/(N*(N-1))
     soma=soma+t
     If(abs(soma-seno)<eps) then
        seno=soma
        exit
     Else
        seno=soma
     End If
             
  End Do
   
End if  
  write(*,*)'Número de termos', N+1

End Function
!!!!!!!!!!!!!!!!!!!!!!
Program testa_Taylor_seno
implicit none
Real:: seno
Real:: x
Real:: f
write(*,*)'Qual o valor de x?'
Read(*,*)x
f= seno(x)
Write(*,*)'Function:',f
Write(*,*)'seno(x):', seno(x)

End program

