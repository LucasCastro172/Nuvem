program cores_imagem
implicit none
integer:: Nx,Ny
Real, allocatable,dimension(:,:):: red,green,blue
Real, allocatable,dimension(:):: f, f1
Integer:: i, n,j,k,s

open(22,file='Red.dat',status='replace',action='write')
open(33,file='Green.dat',status='replace',action='write')
open(44,file='Blue.dat',status='replace',action='write')
open(55,file='fontes.dat',status='old',action='read')
read(55,*)Nx,Ny

write(*,*)Nx,Ny

n=3*Nx*Ny
write(*,*)N

allocate(red(ny,nx),green(ny,Nx),blue(ny,nx))

 Do i=1,nx
  read(55,*)(red(j,i),j=1,ny)
 End Do

 Do i=1,nx
  write(*,*)(red(j,i),j=1,ny)
 End Do

 Do i=1,nx
  read(55,*)(green(j,i),j=1,ny)
 End Do

 Do i=1,nx
  read(55,*)(blue(j,i),j=1,ny)
 End Do
 
 

Do i=1,ny
  write(22,*)(red(i,j),j=1,nx)
  write(33,*)(green(i,j),j=1,nx)  
  write(44,*)(blue(i,j),j=1,nx)
End do

end program
