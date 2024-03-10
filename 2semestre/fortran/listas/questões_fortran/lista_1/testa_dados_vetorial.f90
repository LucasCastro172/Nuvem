program testa_dado 
use jogos
implicit none 
integer:: i,n, Face
integer:: F(6)
real:: t1,t2

call random_seed()
write(*,*)'Quantos lançamentos?'
Read(*,*) N

F=0
call cpu_time(t1)

 Do i=1,N
    call dado_s(face)
    F(face)=F(face)+1
 End do
 
call cpu_time(t2)

write(*,*)'Tempo',t2-t1
 
 Do i=1,6
    !write(*,'("F",i0,":",f10.6)') i,real(F(i))/n
   write(*,*)i,real(F(i))/n
 End do    
     
!Do i=1,N
!   call dado_s(F)
!   if (F==1) then
 !    F1=F1+1
 !  else if (F==2) then
!     F2=F2+1
!   else if (F==3) then
!     F3=F3+1 
!   else if (F==4) then
!     F4=F4+1
!   else if (F==5) then
!     F5=F5+1
!    Else
!     F6=F6+1
!   End if
!End do
!Do i=1,6
! write(*,'("F",i0,":",)')
End program               
