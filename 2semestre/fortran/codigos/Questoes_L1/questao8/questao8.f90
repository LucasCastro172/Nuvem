!programa para testar todas as possibilidades dos lados do triangulo pitagorico
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021
character(3) function trian_ret(a,b,c)
implicit none

integer, intent(in) :: a, b, c
real :: h, c1, c2


	if (a < b + c .and. b  < a + c .and. c < a + b .and. a > 0. .and. b > 0. .and. c > 0) then
		
		if (a > b .and. a > c) then
		h = real(a)
		c1 = real(b)
		c2 = real (c)
		else  if (b > a .and. b > c) then
		h = real(b)
		c1 = real(a)
		c2 = real (c)
		else  if (c > a .and. c > b) then
		h = real(c)
		c1 = real(a)
		c2 = real (b)	
		end if
		
		if (h == sqrt((c1**2) + (c2**2)) )  then
			trian_ret  = 'sim'
		else
			trian_ret = 'nao'
		end if
	
	else
		trian_ret = 'err'
	
	end if
end function


Program Triangulo_ret


Implicit none
Integer :: a, b, c	! laterais do triangulo a, b, c
logical :: cond_1, cond_2, cond_3
integer :: x, i, y, w
integer :: n
character(3), external :: trian_ret
character(3) :: tri 



w= 0
Do x=1,500
	Do y=1,500   
		Do i=1,500 
	     	a = x
	     	b = y
	     	c = i
		tri = trian_ret(a,b,c)
			if( tri == 'sim') then
				if( b > c .and. b > a .and. a > c) then
					write(*,*) a, '', b, '', c
					w = w + 1
				end if
			end if
		end do
	end do
end do

write(*,*) 'o numero de trio pitagoricos é', w



end program
