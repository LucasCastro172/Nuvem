!!!!!!!!!
!Programa para ordenar o maior, o segundo e o terceiro valor colocando a linha correta
!Aluno: Lucas de Castro Costa
!Data: 27/05/2021
program encontra_maior
implicit none
real :: x, maior
real :: y, z
real :: segundo, terceiro
integer :: erro
integer :: i, linha_1
integer :: linha_2, linha_3

open ( 33, file='aa.dat', status='old', action='read')

read(33,*,iostat=erro) x
!read(33,*,iostat=erro) y
!read(33,*,iostat=erro) z
if ( erro == 0 ) then
	segundo = 0.
	terceiro = 0.
	maior = x
	i = 1
	linha_1 = 1
	do
		read(33,*,iostat=erro) x
		if ( erro == 0 ) then
			i = i + 1
			if (x > maior) then
				maior = x
				linha_1 = i
			else if (maior > x .and. x > segundo) then
				segundo = x
				linha_2 = i
			else if (segundo > x .and. x > terceiro) then
				terceiro = x
				linha_3 = i	
			end if
		else
			exit
		end if
	end do
	write(*,*) 'O arquivo tem', i, 'linhas'
	write (*,*) 'O maior valor é', maior, 'na linha', linha_1
	write (*,*) 'O segundo maior valor é', segundo, 'na linha', linha_2
	write (*,*) 'O terceiro maior valor é', terceiro, 'na linha', linha_3
	
else
	write(*,*) 'Erro na primeira leitura!'
end if






end program
