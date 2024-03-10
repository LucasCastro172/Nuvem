!!!!!!!!!
!Programa para calcular o desvio padrão sem utilizar o Rewind
!Aluno: Lucas de Castro Costa
!Data: 25/05/2021
program Media
implicit none
integer :: N
real :: x
integer :: i
real :: m, med, soma, soma_2
real :: desvio, dp

open ( 22, file='dados.dat', status='old', action='read' )

read(22,*) N

soma = 0.
Soma_2 = 0.
do i = 1, N
	read(22,*) x
	soma = soma + x
	soma_2 = soma_2 + x**2
end do

m = soma / N

dp = sqrt( soma_2/n - (soma/N)**2 )

write(*,*) 'Desvio padrão:', dp









end program
