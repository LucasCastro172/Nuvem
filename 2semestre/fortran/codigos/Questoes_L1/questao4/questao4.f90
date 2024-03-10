!Programa para calcular a distancia entre dois pontos
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021
real function distancia (x1,x2,y1,y2)
implicit none
real, intent(in) :: x1, y1 !coordenadas do primeiro ponto
real, intent(in) :: x2, y2 !coordenadas do segundo ponto

distancia = sqrt((x2-x1)**2+(y2-y1)**2) !equacao pra distancia


end function distancia


!Inicio do program

Program Pontos
real, external :: distancia !funcao
real :: x1, y1 !coordenadas do primeiro ponto
real :: x2, y2 !coordenadas do segundo ponto

write(*,*) 'Coordenadas do primeiro ponto'
read(*,*) x1, y1

write(*,*) 'Coordenadas do segundo ponto'
read(*,*) x2, y2

write(*,*) 'Distância entre os dois pontos é', distancia(x1,x2,y1,y2)




end program

