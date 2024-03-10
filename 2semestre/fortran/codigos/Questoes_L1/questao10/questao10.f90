program lancadado
implicit none
integer :: Face
integer :: f1,f2,f3,f4,f5,f6
integer :: i, N



call random_seed()

f1=0
f2=0
f3=0
f4=0
f5=0
f6=0

write(*,*) 'Quantos lançamentos?'
read(*,*) N

do i = 1, N
	Face = dado6()
	if ( Face == 1) then
		f1 = f1 +1
	else if ( Face == 2) then
		f2 = f2 + 1
	else if ( Face == 3) then
		f3 = f3 + 1  
	else if (face == 4) then
		f4 = f4 + 1 
	else if (face == 5) then
		f5 = f5 + 1
	else
		f6 = f6 + 1
	end if
end do

write(*,*) 'F1:', 100 * real(F1) / real(n), '%' 
write(*,*) 'F2:', 100 * real(F2) / real(n), '%'  
write(*,*) 'F3:', 100 * real(F3) / real(n), '%' 
write(*,*) 'F4:', 100 * real(F4) / real(n), '%' 
write(*,*) 'F5:', 100 * real(F5) / real(n), '%' 
write(*,*) 'F6:', 100 * real(F6) / real(n), '%' 




contains

integer function dado6()
implicit none
real :: x

call random_number( x )

dado6 = 6*x + 1.


end function dado6




end program
