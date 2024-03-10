Program jogo_aposta
implicit none
Real:: deposito,aposta,total
character(17):: resultado

deposito=1000.

Do
    Write(*,*)'Valor da aposta:'
    Read(*,*) Aposta

     If(deposito<aposta)then
        Do 
          Write(*,*)'Valor da aposta:'
          Read(*,*) Aposta
           If(deposito>=aposta) then
              Exit
           End if
        End do
     End If

   resultado=craps()
   write(*,*)resultado
  
  If(resultado=='O jogador perdeu.') then
      deposito=deposito-aposta
      write(*,'(A,F8.3)')'depsosito atual:',deposito
  Else
      deposito=deposito+aposta
      write(*,'(A,x,F8.3)')'deposito atual:',deposito
  End if
  
  If(deposito<=0) then
     write(*,*)'Lamento, você faliu'
     Exit
  End if
End Do


contains

integer function dado6()
implicit none
real :: x

call random_number( x )

dado6 = 6*x + 1.

!dado6 = 11*x + 2.


end function dado6


character(17) function craps()
Implicit none
integer:: i,soma,aux

soma=dado6()+dado6()
aux=soma
	If(soma==7.or.soma==11) then
		write(*,*) soma
   		craps='O jogador ganhou!'
	Else if(soma==2.or.soma==3.or.soma==12) then
		write(*,*) soma
         	craps='O jogador perdeu.'
	Else
    		Do
    			write(*,*) soma
    			read(*,*)
       		soma=dado6()+dado6()
        			If(soma==7) then
        				write(*,*) soma
           				craps='O jogador perdeu.'
            			exit
            			end if
        			If(soma==aux) then  
        				write(*,*) soma  
           				craps='O jogador ganhou!'
           			exit  
           			end if                  
    		End Do
	End if       

end function craps

end program 
          
