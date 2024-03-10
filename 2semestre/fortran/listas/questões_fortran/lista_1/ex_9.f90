Real function logax(a,x)
Implicit none
Real, intent(in):: a ! base
Real, intent(in):: x !logaritmando 

If(0<a.and.a/=1.and.x>0)then
   logax=log10(x)/log10(a)
Else 
   write(*,*)'ta errado'
End if  

End function
!===================
Program logaritmo
implicit none
Real:: x, a
Real:: logax

Write(*,*)'Entre com os valores de a,x:'
Read(*,*) a,x

write(*,*) logax(a,x)

End program 

