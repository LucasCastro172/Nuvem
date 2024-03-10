Program teste_U
use Sistemas
implicit none 
integer:: N, g,nc
real(prc), allocatable:: x(:,:), y(:),ys(:),deltay(:)
real(prc), allocatable:: xtx(:,:), xty(:), p(:)
integer:: i

n=100
g=3
nc=g+1 

allocate(x(n,nc),y(n),ys(n),xtx(nc,nc),xty(nc),p(nc),deltay(n))

x(:,1)=1.d0
open(55, file='dados_para_mq.dat', status='old',action='read')
do i=1,n
   read(55,*) x(i,2), y(i)
end do
 close(55)
 
 do i = 3, nc
   x(:,i) = x(:,2)**(i-1)
 end do        
 call Monta_sistema_MQ(N, Nc, x, y, xtx, xty)
 !call fatora_LU(n,xtx)
 !call resolve_LU(N,xtx,xty)
 call elimina(nc, xtx, xty)
 call retrosubs(nc, xtx, xty, p)
 write(*,*) 'A:', p(1)
 write(*,*) 'B:', p(2)
 write(*,*) 'C:', p(3)
 write(*,*) 'D:', p(4)
 ys=matmul(x,p)
open(22,file='P2.out',status='replace')
do i=1,N
  deltay(i) = ys(i) - y(i)
!  write(*,*) x(i,2), y(i), ys(i)
End do
do i=1,N
  write(22,*) x(i,2), y(i), ys(i),deltay(i)
End do

End program 
