program laplace2d
!..... regiao quadrada
implicit none
integer           :: n
real, allocatable :: fi(:,:)
real              :: f1, f2, f3, f4
integer           :: i,j,k
real              :: eps, aux, maior, dif
integer           :: niter

open(44, file='modelo_quadrado.dat',status='old', action='read')
read(44,*) N
read(44,*) f1, f2, f3, f4
read(44,*) eps
 close(44)
 
allocate(fi(n,n))
!preenchendo a matriz
Fi = (f1+f2+f3+f4)/4.
fi(:,1) = f1
fi(1,:) = f2
fi(:,n) = f3
fi(n,:) = f4
!fi(1,1) = (f1+f2)*0.5
!fi(1,n) = (f3+f2)*0.5
!fi(n,n) = (f3+f4)*0.5
!fi(n,1) = (f1+f4)*0.5
niter = 0.
do 
	maior = 0.
	do j = 2, n-1
		do i = 2, n-1
			aux = (fi(i+1,j)+fi(i,j+1)+fi(i-1,j)+fi(i,j-1))*0.25
			dif =abs( aux -fi(i,j))
			if ( dif  > maior ) maior = dif
			fi(i,j) = aux
		end do
	end do
	niter = niter + 1
	if (maior < eps) exit	
end do
write(*,*) 'n inter', Niter

open(55, file='sol_laplace2.dat', status='replace',action='write')

do i = 1,n
	write(55,*) (fi(i,j), j=1,N)
end do


end program
