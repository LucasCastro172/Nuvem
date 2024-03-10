Program jogo_aposta
use jogos
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
End program
