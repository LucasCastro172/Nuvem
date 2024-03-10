Real Function distancia(x1,x2,y1,y2)
Implicit none
Real :: x1,x2,y1,y2
 
distancia=sqrt((x2-x1)**2+(y2-y1)**2)

End Function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program distancia_entre_dois_pontos
Implicit none
Real:: x1,x2,y1,y2,distancia

write(*,*)'Entre com a posição do ponto 1:'
Read(*,*) x1,y1

write(*,*)'Entre com a posição do ponto 1:'
Read(*,*) x2,y2

write(*,*)distancia(x1,x2,y1,y2)

End program
