program laplace2d
!..... regiao quadrada
implicit none
integer           :: n
real, allocatable :: fi(:,:)
real              :: f1, f2, f3, f4
integer           :: i,j,k
integer           :: niter

open(44, file='modelo_quadrado.dat',status='old', action='read')
read(44,*) N
read(44,*) f1, f2, f3, f4
read(44,*) niter
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

do k = 1, niter

	do j = 2, n-1
		do i = 2, n-1
			fi(i,j) = (fi(i+1,j)+fi(i,j+1)+fi(i-1,j)+fi(i,j-1))*0.25
		end do
	end do

end do

open(55, file='sol_laplace.dat', status='replace',action='write')

do i = 1,n
	write(55,*) (fi(i,j), j=1,N)
end do


end program
