program Beni
use estatistica
implicit none
integer:: Nx,Ny, N2x,N2y
Real, allocatable,dimension(:,:):: r,g,b,r2,g2,b2
Real, allocatable,dimension(:):: f, f1
Integer:: i, n,j,k,s

open(22,file='Red.dat',status='replace',action='write')
open(33,file='Green.dat',status='replace',action='write')
open(44,file='Blue.dat',status='replace',action='write')
open(66,file='Red2.dat',status='replace',action='write')
open(77,file='Green2.dat',status='replace',action='write')
open(88,file='Blue2.dat',status='replace',action='write')
open(55,file='Beni.dat',status='old',action='read')
read(55,*)Nx,Ny

write(*,*)Nx,Ny

n=Nx*Ny
write(*,*)N

nx2 = nx/2

ny2 = ny/2

allocate(r(ny,nx),g(ny,Nx),b(ny,nx),r2(ny2,nx2),g2(ny2,Nx2),b2(ny2,nx2))

 Do i=1,nx
  read(55,*)(r(j,i),j=1,ny)
 End Do


! Do i=1,2,nx
!    r2(j,i)=media(4,r(j:j+1,i:i+1))
!    write(66,*)(media(4,r(j:j+1,i:i+1)),j=1,ny)
! End Do


 Do i=1,nx
  read(55,*)(g(j,i),j=1,ny)
 End Do
! Do i=1,2,nx
!   ! r2(j,i)=media(4,r(j:j+1,i:i+1))
!   write(77,*)(media(4,g(j:j+1,i:i+1)),j=1,ny)
! End Do


 Do i=1,nx
  read(55,*)(b(j,i),j=1,ny)
 End Do

do i= 1, nx2
	do j = 1, ny2
		k = 2*i - 1
		s = 2*j - 1
		r2(i,j) = (r(k,1) + r(k+1,s) + )



Do i=1,ny
  write(22,*)(r(i,j),j=1,nx)
  write(33,*)(g(i,j),j=1,nx)  
  write(44,*)(b(i,j),j=1,nx)
End do
end program
