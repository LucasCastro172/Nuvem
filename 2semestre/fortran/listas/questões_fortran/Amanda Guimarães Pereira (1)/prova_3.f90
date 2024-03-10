! Aluna: Amanda Guimarães Pereira
Program dados_circ
use Sistemas
implicit none 
!-----------------------------------------------------------------
! x0: coordenada do eixo x do centro do circulo
! y0: coordenada do eixo y do centro do circulo
! R: raio
! coord_x: coordenada do eixo x
! coord_y: coordenada do eixo y
! p: vetor da solução
! x: matriz dos coeficientes
! y: vetor dos termos indepentes
! xtx: é a matriz resultado da multiplicação da tranposta de x por x.
! xty: é o vetor resultado da multiplicação da tranposta de x por y.
!-------------------------------------------------------------------
real(prc):: x0,y0,R
real(prc), allocatable:: x(:,:), y(:), coord_x(:),coord_y(:)
real(prc), allocatable:: xtx(:,:), xty(:), p(:)
integer:: i,j


allocate(x(15,3),y(15),xtx(3,3),xty(3),p(3),coord_x(15),coord_y(15))

open(55, file='dados_circ.dat', status='old',action='read')
open(66, file='resultado.dat', status='replace')

do i=1,15
   read(55,*) coord_x(i),coord_y(i)
end do
 close(55)

do i=1,15
x(i,1)=coord_x(i) 
x(i,2)=coord_y(i)
x(i,3)= 1.d0
y(i)= coord_x(i)**2 +coord_y(i)**2
end do 

 call Monta_sistema_MQ(15, 3, x, y, xtx, xty)
 call elimina(3, xtx, xty)
 call retrosubs(3, xtx, xty, p)

x0=p(1)/2
y0=p(2)/2
R=sqrt(p(3)+x0**2+y0**2)

write(66,*) 'x0:', x0
write(66,*) 'y0:', y0
write(66,*) 'R:' , R
End program 
