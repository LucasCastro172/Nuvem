program lancadado
implicit none
integer :: Faces
integer :: i, N
integer :: a, b, c, d, e
logical :: cond_v, cond_d
integer :: soma

write(*,*) 'Quantos lançamentos?'
read(*,*) N

soma = dado6() + dado6()
write(*,*) soma
	if ((soma == 7) .or. (soma ==11))  then
	write(*,*) 'Jogador ganhou'
	else if ((soma == 2) .or. (soma == 3) .or. (soma == 12)) then 
	write(*,*) 'Jogador perdeu'
	else if (soma == 4) then
		do i = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 4)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	else if (soma == 5) then
		do a = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 5)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	else if (soma == 6) then
		do b = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 6)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	else if (soma == 8) then
		do c = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 8)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	else if (soma == 9) then
		do d = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 9)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	else 
		do e = 1, N
		soma = dado6() + dado6()
		write(*,*) soma
		if (soma == 10)  then
	        write(*,*) 'Jogador ganhou'
	        stop
	        else if (soma == 7) then 
	        write(*,*) 'Jogador perdeu'
	        stop
	        else 
	        continue
	        end if
		end do
	end if
!write(*,*) 'Vitorias:', f1
!write(*,*) 'Derrotas:' , f2
!write(*,*) 'Continue:' , f3

contains

integer function dado6()
implicit none
real :: x

call random_number( x )

dado6 = 6*x + 1.

!dado6 = 11*x + 2.


end function dado6





end program
