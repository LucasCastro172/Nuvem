! Subrotina utilizada para medir o campo Gravimétrico
! causados por um conjunto de cilíndros

Module cilindro
implicit none
integer, parameter:: prc = kind(1.d0)
real(prc):: a, p, r,  l

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Sub-rotina para integrar
!intervalo de integração (a,b), "e" critério de parada
!soma de todas as contribuições dos trapézio
!Função de entrada
!n = número de trapézios

subroutine trapezios_mod(f, a, b, e, s)              
implicit none
real(prc), intent(in)::a,b, e                         
real(prc), intent(out):: s                            
real(prc), external:: f                               
integer::i, n                                         
real(prc):: h,x
real(prc):: soma, a1, a2,s1, s2, F_ab  
         
         
f_ab = f(a) + f(b)                  
N = 4                                
h = (b-a)/n                        
X = A                               
soma = 0.                          
                                        
do i = 1,n-1                        
	x = x + h                        
	soma = soma + f(x)              
end do                              
                                             
a1= h*( 2*soma+ f_ab )/2.           
         
x = a + h/2  
soma = soma + f(x)               
                                          
do i=1,n-1                       
	x = x + h                   
	soma = soma + f(x)            
end do                          
         
a2 = h*( 2*soma  + f_ab )/4.
     
s1 = (4*a2 - a1)/3.
a1 = a2

do    
	h = h/2.                         
	n = 2*n 
		 
	x = a + h/2    
	soma = soma + f(x)               
                                             
	do i=1,n-1                       
		x = x + h                    
		soma = soma + f(x)            
	end do                           
                                             
	a2 = h*( 2*soma  + f_ab )/4.             
            
    s2 = (4*a2 - a1)/3.
            
     
    if ( abs(s2-s1)> e) then         
		s1 = s2
		a1 = a2
                                 
    else                         
		s = s2                   
		exit
    end if      
end do
         
end subroutine trapezios_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!função Iy
real(prc) function fy ( x )                 
implicit none
real(prc), intent(in):: x
real(prc)::g, h



g = ((p+l)**2+(a-x)**2)*(sqrt((r**2)-(x**2))+sqrt(p**2+r**2+a**2-2*x*a))**2
h = (p**2+(a-x)**2)*(sqrt(r**2- x**2)+sqrt((p+l)**2 + r**2 + a**2 - 2*x*a))**2
!--------------------------------------------------
fy = log(g/h)
	

end function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Monta_sistema_MQ(N, Nc, x, y, xtx, xty)
implicit none
!N =número de pontos, Nc=colunas da matriz x que é um a mais que o grau do polinômio.
integer, intent(in):: N, Nc
real(prc), intent(in):: x(N,Nc), y(N)
real(prc), intent(out):: xtx(Nc,Nc), xty(Nc)
integer:: i, j

xtx = matmul(transpose(x), x)  !Jeito 1
xty = matmul(transpose(x), y)  !menor quantidade de código/pode ser mais rápido

!do i = 1, Nc                         !Jeito 2
!	do j = 1, Nc
!		xtx(i, j) = dot_product(x(:,i), x(:, j))
!		xtx(j, i) = xtx(i, j)
!	end do
!	xty(i) = dot_product(x(:, i), y)
!end do                              ! 
 


end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elimina_pivoteamento(N, A, b, erro)
implicit none 
!Essa subrotina pressupõe que a matriz A seja não singular

integer, intent(in):: N
real(prc), intent(inout):: A(N,N), b(N)
real(prc):: m
integer:: i, j, k, I_max
integer, intent(out):: erro
real(prc)::Maior, aux
real(prc)::Linha(N) 
real(prc):: eps = 1.d-6

erro = 0
do k=1, N
	Maior = A(K, K)
	I_max = k
	do i = k+1, N
		if (abs(A(i, K)) > Maior) then	
			Maior = abs(A(i, k))
			I_max = i
		end if
	end do
	if (Maior > eps) then	
		if (I_max/=k) then
			Linha(k:N) = A(K, K:N)
			A(K, K:N) = A(I_max,K:N)
			A(I_max,K:N) = Linha(K:N)
			aux = b(k)
			b(k) = b(I_max)
			b(I_max) = aux
		end if

		do i = k+1, N
			m = A(i, k) / A(k, k)
			do j = k+1, N                      !esse laço pode ser substituído 
				A(i, j)	= A(i, j) - m*A(k, j)  !substituído pela versão vetorial
			end do                             ! 
			!A(i, k+1:N) = A(i, k+1:N) - n*A(K, K+1:N) !Versão vetorial do laço 
														!gasta mais memória.
			b(i) = b(i) - m*b(k)
		end do 
	else
		erro = 1
		Write(*,*)'Matriz singular ou quase'
		exit
	end if
end do

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine retrosubs(N,U, b, x)
implicit none

!Essa subrotina pressupõe que a matriz U seja triângular superior

integer, intent(in):: N
real(prc), intent(in):: U(N,N)
real(prc), intent(in):: b(N)
real(prc), intent(out)::x(N)
integer:: i, j
real(prc):: soma

x(N) = b(N)/ U(N,N)

do i = N-1, 1, -1
	soma =  0._prc  !0.d0  =  0.prc são iguais
	do j = i+1, N
		soma = soma + U(i,j)*x(j)
	end do
	x(i) = (b(i)-soma)/U(i,i)
end do



end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program cilindros_grav
use cilindro
implicit none
!r = raio do cilíndro
!p = distancia do cilindro até o eixo x
!l = comprimento do cilíndro
!a = ponto onde se quer medir a variação 
real(prc):: cent, eps
integer:: tm 
real(prc), allocatable:: t(:) ,fx(:), y(:), x(:,:), cons(:)
real(prc), allocatable:: xtx(:,:), xty(:), ys(:)
integer:: N, i, j, Nc, erro

eps = 1.e-6 
Nc= 3
N = 40
allocate(t(N), fx(N), y(N), x(N, Nc), cons(Nc))

allocate(xtx(Nc,Nc), xty(Nc), ys(N))

open(11, file='grav_cilindros.dat', status='old', action='read')
do i =1, N
	read(11,*)t(i), y(i)
end do
 close(11)

open(22, file='medidas.dat', status = 'old', action='read')
!read(22,*)
read(22,*)r, l, p
 close(22)

do j = 1, Nc
	write(*,*)'Centro do cilindro'
	read(*,*)cent
	do i = 1,N
		a = t(i) - cent 
		call trapezios_mod(fy, -r, r, eps, fx(i))
		!fo = fx
		
	end do
	
	do i = 1, N
		x(i, j) = fx(i)
	end do
end do	



call Monta_sistema_MQ(N, Nc, x, y, xtx, xty)

call elimina_pivoteamento(Nc, xtx, xty, erro)

call retrosubs(Nc, xtx, xty, cons)

ys = matmul(x, cons)

open(33, file='gra.dat', status = 'replace', action='write')


do i = 1, N
	write(33,*)t(i), y(i), ys(i)
end do

write(*,*) cons(1), cons(2), cons(3)


end program
