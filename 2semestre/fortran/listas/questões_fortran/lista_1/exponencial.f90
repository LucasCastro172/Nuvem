Real function exponencial(x)
! Essa função calcula a função exponencial utilizando a serie de Taylor. 
implicit none 
Real, intent(in):: x
Integer:: N ! número de elementos
Real:: soma
Real:: t  ! termos
Real,parameter:: eps=1.e-4   !precisão para o laço parar 

exponencial=1.
soma=1.
t=1.
n=0
   Do 
      N=N+1
      t=t*x/N  
      soma= soma + t
      If(abs(soma-exponencial)<eps) then
         exponencial= soma
         exit
      else
          exponencial= soma
      End if 
   End Do
 write(*,'(A,x,I2)') 'Números de termos:',N+1
End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program testa_Taylor
implicit none
Real:: exponencial 
Real:: x
Real:: y
write(*,*)'Qual o valor de x?'
Read(*,*)x
y= exponencial(x)
Write(*,*)'Function:',y
Write(*,*)'exp(x):', exp(x)

End program

