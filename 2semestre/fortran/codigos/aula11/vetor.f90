!!!!!!!!!
!Programa para utilizar um vetor
!Aluno: Lucas de Castro Costa
!Data: 01/06/2021
program media_desvio
implicit none
integer :: N
real, allocatable :: v(:)  !vetor
real :: x
integer :: i
real :: m, dp


open ( 22, file='dados.dat', status='old', action='read' )

read(22,*) N

allocate( v(N))

do i = 1, N
	read(22,*) v(i)
end do

 close(22)

call calcula_media_dp(n, v, m, dp)

write(*,*) 'Desvio padrão:', dp

contains

subroutine calcula_media_dp(n, x, m, dp)
implicit none
integer, intent(in) :: N
real, intent(in) :: x(N)
real, intent(out) :: m, dp
real :: soma, soma_2


soma = 0.
Soma_2 = 0.
do i = 1, N
	soma = soma + x(i)
	soma_2 = soma_2 + x(i)**2
end do

m = soma / N

dp = sqrt( soma_2/n - m**2 )


end subroutine









end program
