program teste_ruido
implicit none 
integer, parameter:: N=100
real,dimension(N):: x,y,r
real:: dx
integer:: i,j 

x(1)=-10.
y(1)=sin(x(1))
x(n)=10.
!y(N)=sin(x(n))
y(N)=fun(x(n))
dx=(x(n)-x(1))/(N-1)

Do i=2, N
   x(i) = x(i-1) + dx
!   y(i)=sin(x(i))
    y(i)=fun(x(i))    
End Do

call random_number(r)
r=0.2*(r-0.5)

open(33,file='funcao_ruido.dat',status='replace')
Do i=1,N
  write(33,*) x(i),y(i),y(i)+r(i)
End do   

contains
real function fun(x)
implicit none 
real, intent(in):: x
fun = exp(-abs(x)) * sin(3*x)
end function fun 
          
End program 
