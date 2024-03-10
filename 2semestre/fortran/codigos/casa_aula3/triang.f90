!Dever de casa Aula 3
!Programa para calcular a area de um triangulo com tres lados diferentes
!Aluno: Lucas de Castro Costa
!Data:27/04/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Criar uma funcao para a area do triangulo 
Real function triangulo (a,b,c) 
implicit none
real, intent(in) :: a,b,c !lados do triangulo
!real, intent(in) :: p !metade do perimetro do triangulo
!formula para o calculo da area do triangulo

triangulo = sqrt((a + b + c)/2*((a + b + c)/2 - a)*((a + b + c)/2 - b)*((a + b + c)/2 - c))
!para a realizacao deste trabalho, foi utilizada a formula de Heron

end function triangulo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!iniciar o programa
program areatriangulo 
implicit none
real, external :: triangulo
real :: x,y,z !entrada dos parametros dos lados do triangulo
!real :: p !entrada da metade do perimetro do triangulo
!lados do triangulo
write(*,*) 'Entre com os valores dos lados do triangulo'
read(*,*) x,y,z 

write(*,*) 'area do triangulo e', triangulo(x,y,z)

end program

!metade 
 

