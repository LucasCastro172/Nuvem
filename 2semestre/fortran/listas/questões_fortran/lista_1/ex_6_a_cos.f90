Real function cosseno(x)
implicit none 

Real, intent(in):: x
Integer:: N
Real:: soma
Real:: t
Real, parameter:: eps= 1.e-6 

cosseno=1.
soma=1.
t=1.
N=0.
  
  Do
     N=N+1
     t=-t*x**2/(2*N*(2*N-1))
     soma=soma+t
     If(abs(soma-cosseno)<eps) then
        cosseno=soma
        exit
     Else
        cosseno=soma
     End If
             
  End Do
  
  write(*,*)'Número de termos', N+1

End Function
!!!!!!!!!!!!!!!!!!!!!!
Program testa_Taylor_cosseno
implicit none
Real:: cosseno
Real:: x
Real:: y
write(*,*)'Qual o valor de x?'
Read(*,*)x
y= cosseno(x)
Write(*,*)'Function:',y
Write(*,*)'cosseno(x):', cosseno(x)

End program

