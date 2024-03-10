program numeros_primos
implicit none 
Integer:: i,j, N,primos,k
Real,Allocatable, Dimension(:):: y

write(*,*)'entre com a dimensão'
Read(*,*)n

allocate(y(n))

y=1

primos=0
Do i=2,n
    Do j=i+1,n
       If(y(j)==0.or.mod(j,i)==0) then
         y(j)=0
       Else    
         y(j)=1     
       End if
    End Do  
End Do

!primos=0

y(1)=0.    

Do i=1,n
   If(y(i)==1) then
      primos=primos+1
   End if
End do 
write(*,*)primos
Do i=1,n
write(*,*)y(i)
End Do
    
End program 
