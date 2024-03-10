Program geradados
implicit none 
integer:: N, g
real(8), allocatable:: x(:), y(:), p(:)
real(8):: xi, xf, dx
real(8):: r, ar ! ar: amplitude do ruído. 
integer:: i, j

open(22, file='polinomio_mq.in', status='old', action='read')
Read(22,*)g
Read(22,*)N
Allocate(x(N),y(N),p(g+1))
Read(22,*)p
Read(22,*) xi,xf
Read(22,*)ar

 close(22)
dx=(xf-xi)/(n-1)

do i=1,n
   x(i)= xi + (i-1)*dx
   do j = 1, g+1
      y(i) = y(i) + p(j)*x(i)**(j-1)
   End do
End do 
      
  call random_number(r)
   r =  2*(r - 0.5d0) * ar
   y(i) = y(i) + r

open(33, file='Dados_para_ajustar.dat', status='replace')
  do i=1,n
     write(33,*) x(i), y(i)
  end do
         
End program 
