!!!!!!!!!
!Programa para ordenar
!1. Ler os dados
! - abrir arquivo
! - le o tamanho do vetor
! - alocar a memoria para o vetor
! - le os dados
! 2. ordenar
! 3. saída dos dados
!Aluno: Lucas de Castro Costa
!Data: 03/06/2021
program dados_ordenados
implicit none
integer :: N
real, allocatable :: v(:)  !vetor
integer :: i
real :: md
real :: t1, t2

open ( 33, file='dados.dat', status='old', action='read' )
read(33,*) N
allocate( v(n)) 
do i = 1, n
	read(33,*) v(i)
end do
 close(33)
 
call cpu_time(t1)
call ordenar_bolha(n,v)
call cpu_time(t2)
write(*,*) 'Tempo:', t2 - t1
 
open ( 44, file='ordem.dat', status='replace', action='write' )
do i = 1, N
 write(44,*) v(i)	
end do
 close(44)

contains

subroutine Ordenar_bolha(n , v)
implicit none
! . . .
integer, intent(in) :: N
real, intent(inout) :: v (n)
real :: aux
integer :: i, j, troca

i = 0
do 
	i = i + 1
	troca = 0
	do j = 1, N-1
		if ( v(j) > v(j+1) ) then
			aux = v(j)
			v(j) = v(j+1) 
			v(j+1) = aux
			troca = 1
		end if
	end do
        if (troca == 0) exit
end do


end subroutine ordenar_bolha



end program
