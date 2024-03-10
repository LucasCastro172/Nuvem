!!!!!!!!!!
!Programa para gerar os valores de periodo, resistividade aparente e fase em graus de uma
!sondagem magnetotelurica
!Aluno: Lucas de Castro Costa
!Data: 23/07/2021
program magnetotelurica
implicit none
Integer,parameter   :: pr = kind(1.d0) 
reaL(pr),parameter  :: pi = 3.14159265359
real(pr)            :: per   !permeabilidade
real(pr)            :: eps 
real(pr)            :: rho    !resistividade aparente
real(pr)            :: w !frequencia angular
integer             :: nf !numero de frequencias
real(pr)            :: f,T !frequencia e periodo
complex(pr)         :: z    !impedancias
real(pr)            :: phi, phi_graus !fase(em radianos) 
integer             ::  i
real(pr)            :: x,y !parte real e parte imaginária

eps = 1.d-7
per = 4*pi*eps !permeabilidade

open(11,file='mt1d.out',status='old',action='read') 
open(22,file='sond.out',status='replace',action='write')

read(11,*) nf 

Do i=1,nf 
   read(11,*) f,z  !lendo os valores de frequencia e impedancia
   T = 1./f
   x = real(z) !parte real complexo
   y = aimag(z)! parte imag complexo
   w = 2.*pi*f
   rho = (1./(w*per))*(x**2+y**2) !a raiz representa o absoluto do numero complexo
   phi = atan2(y,x) !calculo do phi(resposta em rad)
   phi_graus = phi * (180/pi) !conversão do phi de rad para gruas
   write(22,*) T, rho, phi_graus !guardando os valores em um arquivo de saída
End do


end program
