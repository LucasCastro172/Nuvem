!!!!!!!!!!
!Programa para calcular a media de um conjunto de numeros
!Aluno: Lucas de Castro Costa
!Data: 25/05/2021
program Media
implicit none
integer :: N
real :: x, soma
integer :: i
real :: m

open ( 22, file='dados.dat', status='old', action='read' )

read(22,*) N

soma = 0.
do i = 1,N
	read(22,*) x
	soma = soma + x
end do

m = soma / N

write(*,*) 'Média', m



end program
