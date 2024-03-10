Module jogos
implicit none

Contains
 
Integer function dado()
implicit none
Real::x 

call random_number(x)
dado=6*x+1

End function
!!!!!!!!!!!!!!!!!!!!!!!
subroutine dado_s(face)
implicit none
integer, intent(out):: face
Real::x 

call random_number(x)
face=6*x+1

End subroutine dado_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!character(17) function craps()
!Implicit none
!integer:: i,soma,dado1,dado2,aux

!dado1=dado()
!dado2=dado()
!soma=dado1+dado2
!aux=soma
!If(soma==7.or.soma==11) then
!   craps='O jogador ganhou!'
!Else if(soma==2.or.soma==3.or.soma==12) then
!         craps='O jogador perdeu.'
!Else
!    Do
!       dado1=dado()
!       dado2=dado()
!       soma=dado1+dado2
!        If(soma==7) then
!           craps='O jogador perdeu.'
!            exit
!        If(soma==aux) then
!            craps='O jogador ganhou!'
!           exit                    
!        End If
!    End Do
!End if        
                 
!Else
!    Do i=1,7
!       dado1=dado()
!       dado2=dado()
!       soma=dado1+dado2
!       aux=soma
!        If(soma==7) then
!           craps='O jogador perdeu.'
!            exit
!        Else
!            craps='O jogador ganhou!'        
!        End If
!    End Do
!End if        

!End function
End module jogos
