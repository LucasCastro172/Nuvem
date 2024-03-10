!!!!!!!!!
program Produto_Matriz_vetor
implicit none
integer  :: m, n
real :: A(4,3), x(3), y(4)
integer :: i


m = 3
n = 3
A(1,:) = [1. , 2., 3.]
a(2,:) = [2., 1., 2.]
a(3,:) = [1., 2., 1.]
x  = 1.

open(55, file='test1.bin', status='replace',action='write')

do i = 1, m
	write(55,*) A(i,:)
end do

open(66, file='test2.bin', status='replace',action='write')

do i = 1, m
	write(66,*) A(i,:)
end do


!write(*,*) 'Vetor x:'
!write(*,*) x

!call Prod_mat_vet( M, n, A, x , y)
!write(*,*) "Vetor y = Ax"
!write(*,*) y




!contains

!subroutine Prod_mat_vet( M, n, A, x , y)
!!....
!implicit none
!integer, intent(in) :: M !Numero de linhas da matriz
!integer, intent(in) :: N ! Numero de colunas...
!real, intent (in) :: A(m,n), x(n)
!real, intent (out) :: y(m)
!integer :: i, j

!!do i = 1, m
!!	y(i) = dot_product( A(i,:), x )
!!end do

!y=  0
!do i = 1, M
!	do j  = 1, n
!		y(i) = y(i) + A(i,j)*x(j)
!	end do
!end do



!end subroutine Prod_mat_vet



end program
