program raizes
implicit none
Real:: A,B,C ! coeficientes
Real:: s1,s2 ! soluções
Integer:: N  ! número de raízes

write(*,*)
write(*,*)'Cálculo das raízes da equação Ax^2+Bx+C=0'
write(*,*)
write(*,*)'Escreva os valores de A, B, C, nesta ordem'
write(*,*) A, B, C

call eq_2o_gr(A,B,C,s1,s2,N)

If(N==2) then
   write(*,*)
   write(*,*)'Duas raízes reais:'
   write(*,*)'R1:',s1
   write(*,*)'R2:',s2
Else if(N==1) then
   write(*,*)
   write(*,*)'Uma raíz de multiplicidade 2' 
   write(*,*)'R1=R2=',s1
Else if(N==0) then
    write(*,*)   
    write(*,*)'Delta menor que zero. Não há raízes reais'
End if
write(*,*)

End program    
!!!!!!!!!!!!!!!!!!
Subroutine eq_2o_gr(A,B,C,r1,r2,N)
Implicit none
Real,intent(in):: A,B,C
Real,intent(out):: r1,r2
Integer,intent(out):: N 
Real:: Delta

Delta=B**2-4*A*C

If(Delta>0.) then
  r1=(-B - sqrt(Delta))/(2*A)
  r2=(-B + sqrt(Delta))/(2*A)
  N=2
If(Delta<0.) then
  N=0
Else
   r1=-B/(2*A)
   r2=r1
   N=1
End if
   
End Subroutine eq_2o_gr
