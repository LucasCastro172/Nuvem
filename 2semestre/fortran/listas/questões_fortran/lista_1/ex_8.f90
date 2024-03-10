Program triangulos_retangulos
implicit none
Integer:: i,j,k ! contadores
Integer:: l1,l2,l3 ! lados do triângulo
character(32):: resultado,triangulo_retangulo

Open(33,file='triangulos_retângulos.dat',status='replace',action='write')

Do k=1,500
  Do j=1,500   
     Do i=1,500 
      l1=i
      l2=j
      l3=k
      resultado=triangulo_retangulo(l1,l2,l3)
      If(l2>l3)then
         if(resultado=='É um triângulo retângulo') then
             write(33,*)l1,l2,l3
             !write(*,*)resultado
         End if
      End if
     End Do 
  End Do
End Do  
End program     
