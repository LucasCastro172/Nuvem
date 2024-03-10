!!!!!!!!!
!Programa para calcular o desvio padrão sem utilizar o Rewind
!Aluno: Lucas de Castro Costa
!Data: 25/05/2021
program encontra_maior
implicit none
real :: x, maior
real :: segundo, terceiro
integer :: erro
integer :: i, linha

open ( 33, file='aa.dat', status='old', action='read')

read(33,*,iostat=erro) x
if ( erro == 0 ) then
	maior = x
	i = 1
	linha = 1
	do
		read(33,*,iostat=erro) x
		if ( erro == 0 ) then
			i = i + 1
			if (x > maior) then
				maior = x
				linha = i
			end if
		else
			exit
		end if
	end do
	write(*,*) 'O arquivo tem', i, 'linhas'
	write (*,*) 'O maior valor é', maior, 'na linha', linha
else
	write(*,*) 'Erro na primeira leitura!'
end if






end program
