PROGRAM desvio_padrao
implicit none
Integer:: n,i ! número de elementos
Real:: dp,soma,x,d2

Open(33,file='elementos.dat',status='old',action='read')
Read(33,*) n
write(*,*) n

d2=0.
soma=0.
Do i=1,n
  Read(33,*) x
  d2 = d2 + x**2
  soma = soma + x
End do
write(*,*) d2,soma
  ! DO i=1,n
   !  d2(1)=d(1)**2
   !  d2(i+1)=d2(i)+d(i+1)**2
   !  soma(1)=d(1)
   !  soma(i+1)=soma(i)+d(i+1)
   !  write(*,*) d2(i),soma(i)
   !END DO
   dp=sqrt(d2/n-soma**2/n**2)
   write(*,*) dp
END PROGRAM

