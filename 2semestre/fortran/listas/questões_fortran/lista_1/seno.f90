Real function seno(x)
implicit none 

Real, intent(in):: x
Integer:: N
Real:: soma
Real:: t
Real, parameter:: eps= 1.e-6 

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
  
  write(*,*)'Número de termos', N+1

End Function
!!!!!!!!!!!!!!!!!!!!!!
Program testa_Taylor_seno
implicit none
Real:: seno
Real:: x
Real:: y
write(*,*)'Qual o valor de x?'
Read(*,*)x
y= seno(x)
Write(*,*)'Function:',y
Write(*,*)'seno(x):', seno(x)

End program

  
