program Beni
!use estatistica
implicit none
integer:: Nx,Ny, N2x,N2y
Real, allocatable,dimension(:,:):: r,g,b,r2,g2,b2
Integer:: i, n,j,k,s

open(22,file='Red.dat',status='replace',action='write')
open(33,file='Green.dat',status='replace',action='write')
open(44,file='Blue.dat',status='replace',action='write')
open(66,file='Red2.dat',status='replace',action='write')
open(77,file='Green2.dat',status='replace',action='write')
open(88,file='Blue2.dat',status='replace',action='write')
open(55,file='teste.dat',status='old',action='read')
read(55,*)Nx,Ny

write(*,*)Nx,Ny

n=Nx*Ny
write(*,*)N

allocate(r(ny,nx),g(ny,Nx),b(ny,nx),r2(ny/2,nx/2),g2(ny/2,Nx/2),b2(ny/2,nx/2))

 Do i=1,nx
  read(55,*)(r(j,i),j=1,ny)
 End Do
   
   Do i=1,nx
  write(22,*)(r(j,i),j=1,ny)
 End Do


       Do j=2,7,3
            Do i=1,1
              r2(i,j-1)=sum(r(j:j+1,i:i+1)) 
  !r2(1,3)=r(1,3)+r(1,4)+r(2,3)+r(2,4)
!write(*,*)r2(1,1),r2(1,3) 
          End do 
       End do  
  
Do i=1,1
  write(66,*)(r2(i,j),j=1,4)
 ! write(33,*)(g(i,j),j=1,nx)  
  !write(44,*)(b(i,j),j=1,nx)
End do
end program
