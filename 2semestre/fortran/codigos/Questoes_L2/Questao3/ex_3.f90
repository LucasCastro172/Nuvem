program media_movel
use estatistica
!implicit none
integer:: N,w
Real, allocatable, Dimension(:):: d, x
Real, allocatable, Dimension(:):: m, y,v
integer:: i,k

open(22,file='dados_ruido.dat',status='old',action='read')
open(33,file='dados_suavizado_n6.dat',status='replace',action='write')
Read(22,*) N

write(*,*)'Tamanho da janela,w?'
Read(*,*) w

Allocate(x(N),d(N),m(n-w+1),y(n-w+1),v(w))

Do i=1,N
Read(22,*) x(i),d(i)
End do

k=0
Do i=1, N-w+1
m(i)=media(w,d(i:i+w-1))
y(i)=(x(i+3)-x(i+2))/2 + x(i+2)
write(33,*) y(i),m(i)
k=k+1
End do

write(*,*)'Número de laços:'
write(*,*)k

End program 
