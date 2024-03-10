Program teste_U
use Sistemas
implicit none 
integer:: N, g,nc
real(prc):: Y0,fi
real(prc), allocatable:: x(:,:), y(:),ys(:), t(:)
real(prc), allocatable:: xtx(:,:), xty(:), p(:)
integer:: i,j

n=50
g=1
nc=g+1 

allocate(x(n,nc),y(n),ys(n),xtx(nc,nc),xty(nc),p(nc),t(n))

!x(:,1)=1.d0
open(55, file='dados_seno.dat', status='old',action='read')
do i=1,n
   read(55,*) t(i),y(i)
end do
 close(55)

Do i=1,n
x(i,1)=sin(t(i)) 
x(i,2)=cos(t(i))
End do 

Do i=1,n
write(*,*)(x(i,j),j=1,2)
End do 

! do i = 3, nc
!   x(:,i) = x(:,2)**(i-1)
! end do        
 call Monta_sistema_MQ(N, Nc, x, y, xtx, xty)
 !call fatora_LU(n,xtx)
 !call resolve_LU(N,xtx,xty)
 call elimina(nc, xtx, xty)
 call retrosubs(nc, xtx, xty, p)
 write(*,*) 'p1:', p(1)
 write(*,*) 'p2:', p(2)
Y0=sqrt(p(1)**2+p(2)**2)
fi=atan(p(2)/p(1))
 write(*,*) 'Y0:', Y0
 write(*,*) 'fi:', fi
 ys=matmul(x,p)
open(22,file='ex_3.out',status='replace')
do i=1,N
  !write(22,*) x(i,2), y(i), ys(i),deltay(i)
  write(22,*) t(i),y(i), ys(i)
End do

End program 
