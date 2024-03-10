
Program Triangulo_ret

!Programa para determinar se os três valores podem representar um triangulo
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021

Implicit none
Integer :: a, b, c	! laterais do triangulo a, b, c
logical :: cond_1, cond_2, cond_3
	

write(*,*) 'Valores para um possível triângulo retangulo'
read(*,*) a, b, c


cond_1 = (a > 0.) .AND. (b > 0.) .AND. (c > 0.0)
cond_2 = (a + b > c) .AND. (a + c > b) .AND. (b + c > a)
cond_3 = (a**2 == b**2 + c**2) .OR. (b**2 == a**2 + c**2) .OR. (c**2 == b**2 + a**2) 

	If (cond_1 .AND. cond_2 .AND. cond_3) THEN
		
		Write(*,*)  'Estes valores podem representar um triângulo retângulo'
	ELSE
		Write(*,*)  'Estes valores não podem representar um triângulo retângulo'
	END IF

end program Triangulo_ret
