!Só irá existir um triângulo se, somente se, os seus lados obedeceram à seguinte regra: um de seus lados deve ser maior que o valor absoluto (módulo) da diferença 
!dos outros dois lados e menor que a soma dos outros dois lados. Veja o resumo da regra abaixo: 
character(21) function lados_do_triangulo(A,B,C)
Real, intent(in):: A,B,C ! lados
 
 If(abs(A-B)<C.and.C<(A+B).and.abs(A-C)<B.and.C<(A+C).and.abs(B-C)<A.and.A<(B+C)) then
    lados_do_triangulo='É um triângulo'
    !write(*,*)'É um triângulo'
 Else 
    lados_do_triangulo='Não É um triângulo'
    !write(*,*)'Não É um triângulo'
 End if

End function 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

Program triangulo
Real:: A,B,C
character(21)::lados_do_triangulo

Write(*,*)'Entre com o valor dos lados'
Read(*,*) A,B,C

write(*,*) lados_do_triangulo(A,B,C)

End program     
