program gera_dados
implicit none 
Real:: x
Integer:: i,n

call random_seed()

write(*,*) 'Quantos valores?'
read(*,*) N

open(33,file='vetor.dat',status='replace')

write(33,*)N

Do i=1,N
   call random_number(x)
   x=2*(x-0.5)
   write(33,*) x
End Do    

End program
