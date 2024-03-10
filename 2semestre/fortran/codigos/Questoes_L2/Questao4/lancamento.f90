Program lancamento_de_dados
implicit none
Integer:: N,l,i,face
Integer, dimension(:), allocatable:: F

open(33,file='lancamentos.dat',status='replace',action='write')

write(*,*)'Quantidades de dados?'
read(*,*)N

write(*,*)'Quantos lançamentos?'
read(*,*)l

allocate(f(6*n))

call random_seed()
F=0

Do i=1,l
 Call lancamento_dados(n,face)
 F(face)=F(face)+1
End do

Do i=1,6*n
write(33,*) i, real(f(i))/l  
End do 

End program 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine lancamento_dados(n,soma)
implicit none
Integer, intent(in)  :: n
Integer,intent(out)  :: soma
real                 :: x
Integer              :: face(n),i
 
soma=0 
 Do i=1,n  
   call random_number(x)
  face(i)=6*x+1
  soma=soma+face(i)
 End Do
!write(*,*)soma
End subroutine
