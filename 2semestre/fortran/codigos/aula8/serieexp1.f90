program exponencial
implicit none
integer :: N


write(*,*) 'N?'
read(*,*) N

write(*,*) 'Fatorial:', fatorial(N)

Contains

integer function fatorial (N)
implicit none
integer, intent(in) :: N
integer :: i

fatorial = 1
do i = 2, N
	fatorial = fatorial * i
end do

end function fatorial



end program
