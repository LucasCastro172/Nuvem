Program craps
use jogos
Implicit none
integer:: i,soma,dado1,dado2,aux

dado1=dado()
dado2=dado()
soma=dado1+dado2
aux=soma
write(*,*)dado1,dado2,soma

If(soma==7.or.soma==11) then
   write(*,*)'O jogador ganhou!'
Else if(soma==2.or.soma==3.or.soma==12) then
         write(*,*)'O jogador perdeu.'
Else
    Do
       dado1=dado()
       dado2=dado()
       soma=dado1+dado2
        If(soma==7) then
           craps='O jogador perdeu.'
            exit
        If(soma==aux) then
            craps='O jogador ganhou!'
           exit                    
    End Do
End if        

!Else
!    Do i=1,7
!       dado1=dado()
!       dado2=dado()
!       soma=dado1+dado2
!       write(*,*)dado1,dado2,soma
!        If(soma==7) then
!           write(*,*)'O jogador perdeu'
!            exit
!        Else
!            write(*,*)'O jogador ganhou!'
!        End If
!    End Do
!End if        
End program 
