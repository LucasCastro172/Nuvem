program Dados_organizados
implicit none
integer :: N,i,erro,Nx
real, allocatable :: x(:), y(:)


open ( 22, file='coordenadas.dat', status='old', action='read',iostat=erro )

do

	read(22,*,iostat=erro) N
	if ( erro == 0) then
		cycle
	else
		exit
	end if
	
end do

allocate(x(N),y(N))
rewind(22)

do i = 1, n
	read(22,*) Nx, x(i), y(i)
end do
 close(22)
 
call ordenar2d(x,y,n)

open ( 33, file='coordenadas_xy.dat', status='old', action='write')

do i = 1,n
	write(2,*) i, x(i), y(i)
enddo

 close(33)
 
write(*,*) "Programa concluido"


contains

subroutine ordenar2d(x,y,N)
implicit none
integer, intent(in) :: n
real, intent(inout) :: x(N), y(n)
real                :: aux, aux2
integer             :: i,j, troca

i = 0
do 
	i = i  + 1
	troca = 0
	do j = 1, n-i
	
		if ( x(j) > x(j+1)) then
			troca = 1
			aux = x(j)
			aux2 = y(j)
			x(j) = x(j+1)
			y(j) = y(j+1)
			x(j+1)  = aux
			y(j+1) = aux2
		
		end if
		
	end do
	if ( troca == 0) exit
	
end do


end subroutine ordenar2d


end program
