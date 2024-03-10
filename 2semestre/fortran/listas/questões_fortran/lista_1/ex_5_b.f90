character(32) function triangulo_retangulo(A,B,C)
Integer, intent(in):: A,B,C ! lados
 
 If(B<A.and.C<A.and.A**2==(B**2 + C**2)) then
    triangulo_retangulo='É um triângulo retângulo'
    !write(*,*)'É um triângulo retângulo'    
 !Else If(A<B.and.C<B.and.B**2==A**2 + C**2) then
  !  triangulo_retangulo='É um triângulo retângulo'   
 !Else If(A<C.and.B<C.and.C**2== B**2 + C**2) then
 !   triangulo_retangulo='É um triângulo retângulo'    
 Else
    triangulo_retangulo='Não é um triângulo retângulo'
    !write(*,*)'Não é um triângulo retângulo'
 End if

End function 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!Program condicao_triangulo_retangulo
!Integer:: lado_1,lado_2,lado_3
!character(32)::triangulo_retangulo

!Write(*,*)'Entre com o valor dos lados começando pelo maior valor:'
!Read(*,*) lado_1,lado_2,lado_3

!write(*,*) triangulo_retangulo(lado_1,lado_2,lado_3)

!End program     
