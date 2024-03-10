program lancadado
implicit none
integer :: Face
integer :: F(6)
integer :: i, N



call random_seed()

F = 0

write(*,*) 'Quantos lançamentos?'
read(*,*) N

do i = 1, N
	face = dado6()
	f(face) = f(face) + 1
end do

do i = 1,6
	write(*,*) 'Face',i, 100 * real(F(i))/real(n),'%'
end do

!Face = dado6()
!write(*,*) Face
!Face = dado6()
!write(*,*) Face
!Face = dado6()
!write(*,*) Face
!Face = dado6()
!write(*,*) Face
!Face = dado6()
!write(*,*) Face
!Face = dado6()
!write(*,*) Face





contains

integer function dado6()
implicit none
real :: x

call random_number( x )

dado6 = 6*x + 1.


end function dado6




end program
