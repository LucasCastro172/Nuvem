program Produto_Matriz_Matriz
implicit none
integer  :: ma, nab, mb
real :: A(4,4), B(4,3), C(4,3)
integer :: i


ma = 4
nab = 4
mb = 3
!Primeira matriz
A(1,:) = [1., 0., 1., 0.]
a(2,:) = [2., 1., 0., 1.]
a(3,:) = [3., 0., 2., 1.]
a(4,:) = [4., 1., 1., 0.]

!Segunda matriz
B(1,:) = [1., 2., 1.]
B(2,:) = [1., 1., 0.]
B(3,:) = [1., 2., 0.]
B(4,:) = [1., 1., 1.]


write(*,*) 'Matriz A:'
do i = 1, ma
	write(*,*) A(i,:)
end do
write(*,*) 'Matriz B:'
do i = 1, nab
	write(*,*) B(i,:)
end do

call Prod_mat_mat( Ma, Nab, Mb, A, B , C)
write(*,*) "Matriz C = AB"
do i = 1, ma
	write(*,*) C (i,:)
end do




contains

subroutine Prod_mat_mat( Ma, Nab, Mb, A, B , C)
!....
implicit none
integer, intent(in) :: Ma !Numero de linhas da matriz  A
integer, intent(in) :: Nab ! Numero de colunas A e linhas B
integer, intent(in) :: Mb !numero de colunas B
real, intent (in) :: A(ma,nab), B(nab,mb)
real, intent (out) :: C(ma,mb)
integer :: i, j, z


! C = 0
!do i = 1, ma 
!	do j = 1, mb
!	C(i,j) = dot_product( A(i,:), B(:,j) )
!	end do
!end do


 C = 0 
do i  = 1, Ma
	do j = 1, mb
		do z = 1, nab
			C(i,j) = C(i,j) + A(i,z)*B(z,j)
		end do
	end do
end do


end subroutine Prod_mat_mat



end program
