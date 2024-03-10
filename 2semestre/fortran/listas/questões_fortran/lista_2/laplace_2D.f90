Program Eq_laplace_2D
implicit none 
integer:: nx,ny
real:: front(4)
real,allocatable:: Fi(:,:)
real::tol
integer:: i

Open(22,file='laplace_retangulo.in', status='old',action='read')
Read(22,*) nx,ny
allocate(Fi(ny,nx))
read(22,*)Front
read(22,*)tol
 close(22)
 
 call laplace(nx,ny,front,fi,tol)  

open(33,file='laplace.dat',status='replace')
  Do i=1,Ny
    write(33,*)Fi(i,:)
  End Do  

Contains
subroutine laplace(nx,ny,front,fi,eps)
implicit none 
integer, intent(in):: nx,ny
real, intent(in):: front(4)
real, intent(in):: eps
real, intent(inout):: Fi(ny,nx)
real:: dif,dif_max, Fi_novo
integer:: i,j,niter

Fi=sum(front)/4.

Fi(:,1)=Front(1)
Fi(1,:)=Front(2)
Fi(:,Nx)=Front(3)
Fi(Ny,:)=Front(4)

Fi(1,1)=(Front(1)+front(2))*0.5
Fi(1,Nx)=(Front(2)+front(3))*0.5
Fi(Ny,Nx)=(Front(3)+front(4))*0.5
Fi(Ny,1)=(Front(1)+front(4))*0.5

Niter=0

Do 
  dif_max=0.
  do i=2,Ny-1
     do j=2,Nx-1
	Fi_novo=0.25*(Fi(i,j-1))+(Fi(i,j+1))+(Fi(i-1,j))+(Fi(i+1,j))
	dif=abs(Fi_novo-Fi(i,j))
	Fi(i,j)=Fi_novo
	if(dif > dif_max) dif_max=dif
     End do
   End do 
   Niter=Niter+1
   if(dif_max < eps) exit
End do     	
write(*,*)'Niter:', Niter
End subroutine laplace
End program
