!!!!!!!!!!
!Programa para calcular o desvio padrão 
!Aluno: Lucas de Castro Costa
!Data: 25/05/2021
program Media
implicit none
integer :: N
real :: x
integer :: i
real :: m, med, soma
real :: desvio, dp

open ( 22, file='dados.dat', status='old', action='read' )

read(22,*) N

m = 0.
do i = 1,N
	read(22,*) x
	m = m + x
end do

med = m / N

!write(*,*) 'Média', m

rewind(22)

read(22,*) N

desvio = 0.
soma = 0.
do i = 1,N
	read(22,*) x
	desvio = (x - med)**2
	soma = soma + desvio
end do


dp = sqrt(soma/N)


write(*,*) 'Desvio padrão:', dp






end program
