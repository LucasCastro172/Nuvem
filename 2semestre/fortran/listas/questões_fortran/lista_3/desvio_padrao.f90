PROGRAM desvio_padrao
Implicit none
Integer:: n ! número de elementos
Real,Allocatable, Dimension(:)::d,d2,soma ! desvio padrão, elemento,elemento ao quadrado
Real:: dp
Open(44,file='desvio_padrão.dat',status='replace',action='write')
Open(33,file='elementos.dat',status='old',action='read')

Read(33,*) n
write(*,*) n
Allocate(d(n),d2(n),soma(n))

Do i=1,n
  Read(33,*) d(i)
 ! write(*,*) d(i)
End do

   DO i=1,n
     d2(1)=d(1)**2
     d2(i+1)=d2(i)+d(i+1)**2
     soma(1)=d(1)
     soma(i+1)=soma(i)+d(i+1)
     write(*,*) d2(i),soma(i)
   END DO
   dp=sqrt(d2(n)/n-soma(n)**2/n**2)
   write(*,*) dp
END PROGRAM


