program testando
implicit none
integer:: lixo 
Real:: x
Integer:: i
Open(33,file='testando.dat',status='old',action='read')
Read(33,*)
Do i=1,7
Read(33,'(x,F8.3)') x
write(*,*) x
End Do
End program
