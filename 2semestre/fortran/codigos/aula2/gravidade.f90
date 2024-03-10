!segundo codigo
program gravidade
!Programa para calcular a aceleracao da gravidade
!Aluno: Lucas de Castro Costa
!Data:22/04/2021
implicit none
real,parameter :: q = 6.6743e-11 !constante gravitacional  Unidade=Newton metro quadrado por  kilograma ao quadrado
real :: g !aceleracao gravitacional (saida)
real :: m !massa do planeta (entrada)
real :: r !raio do planeta (entrada)

!!!!Entradas
write(*,*) 'Massa do planeta (kg):'
read (*,*) m

write(*,*) 'Raio do planeta (m):'
read(*,*) r

!!!!!!!!!Possibilidade de entrada do raio do planeta em km
!write(*,*) 'Raio do planeta (m):'
!read(*,*) r
!r = r * 1000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!
!equacao da aceleracao gravitacional
g = q*m/r**2
!!!saida
write(*,*) 'A aceleracao gravitacional e (m/s^2)', g
!!!


end program
