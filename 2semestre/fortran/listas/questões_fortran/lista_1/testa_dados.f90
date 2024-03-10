program testa_dado 
use jogos
implicit none 
integer:: i,n, F
integer:: F1,F2,F3,F4,F5,F6
integer:: F(6)
call random_seed()

F=[]
!F1=0
!F2=0
!F3=0
!F4=0
!F5=0
!F6=0

write(*,*)'Quantos lançamentos?'
Read(*,*) N

Do i=1,N
   call dado_s(F)
   if (F==1) then
     F1=F1+1
   else if (F==2) then
     F2=F2+1
   else if (F==3) then
     F3=F3+1 
   else if (F==4) then
     F4=F4+1
   else if (F==5) then
     F5=F5+1
    Else
     F6=F6+1
   End if
End do
write(*,*)'F1:', real(F1)/real(N)
write(*,*)'F2:', real(F2)/real(N)
write(*,*)'F3:', real(F3)/real(N)
write(*,*)'F4:', real(F4)/real(N)
write(*,*)'F5:', real(F5)/real(N)
write(*,*)'F6:', real(F6)/real(N)

End program               
