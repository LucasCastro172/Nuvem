Program Triangulo

!Programa para determinar se os três valores podem representar um triangulo
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021

Implicit none
Real:: a, b, c	! laterais do triangulo a, b, c
logical :: cond_1, cond_2
	

write(*,*) 'Valores para um possível triângulo'
read(*,*) a, b, c


cond_1 = (a > 0.) .AND. (b > 0.) .AND. (c > 0.0)
cond_2 = (a + b > c) .AND. (a + c > b) .AND. (b + c > a)

	If (cond_1 .AND. cond_2) THEN
		
		Write(*,*)  'Estes valores podem representar um triângûlo'
	ELSE
		Write(*,*)  'Estes valores não podem representar um triângulo'
	END IF

end program Triangulo

