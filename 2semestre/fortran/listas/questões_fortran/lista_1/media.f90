program media
implicit none
real:: x,soma, md
integer:: i,n

open(22,file='elementos.dat',status='old',action='read')
read(22,*)n

soma=0.
Do i=1,n
  Read(22,*)x
  soma=soma+x
End do 

md=soma/n

Rewind(22)

write(*,*)'media',md
write(*,*)'desvio padrão',dp

End program  
