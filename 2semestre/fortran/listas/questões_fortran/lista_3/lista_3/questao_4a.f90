module Grav_cilindro
use Sistemas
implicit none
Real(prc) :: R, p, L, a
!integer, parameter:: prc = kind(1.d0)

Contains
!-----------------------------------
Real(prc) function funcao(x)
implicit none
Real(prc), intent(in):: x
Real(prc) :: nu, d

nu =((p+L)**2+(a-x)**2)*(sqrt(R**2-x**2)+sqrt(p**2+R**2+ a**2-2*x*a))**2
d =(p**2+(a-x)**2)*(sqrt(R**2-x**2)+sqrt((p+L)**2+R**2+ a**2-2*x*a))**2
funcao=log(nu/d)

End Function
!-----------------------------------
subroutine trapezios(f, a, b, e, s)              
implicit none
real(prc), intent(in)::a,b, e                         
real(prc), intent(out):: s                            
real(prc), external:: f                               
integer::i, n                                         
real(prc):: h,x
real(prc):: soma, s1n, s2n,s1, s2, F_ab  
         
         
f_ab = f(a) + f(b)                  
N = 4                                
h = (b-a)/n                        
X = A                               
soma = 0.                          
                                        
do i = 1,n-1                        
	x = x + h                        
	soma = soma + f(x)              
end do                              
                                             
s1n= h*( 2*soma+ f_ab )/2.           
         
x = a + h/2  
soma = soma + f(x)               
                                          
do i=1,n-1                       
	x = x + h                   
	soma = soma + f(x)            
end do                          
         
s2n = h*( 2*soma  + f_ab )/4.
     
s1 = (4*s2n - s1n)/3.
s1n = s2n

do    
	h = h/2.                         
	n = 2*n 
		 
	x = a + h/2    
	soma = soma + f(x)               
                                             
	do i=1,n-1                       
		x = x + h                    
		soma = soma + f(x)            
	end do                           
                                             
	s2n = h*( 2*soma  + f_ab )/4.             
            
    s2 = (4*s2n - s1n)/3.
            
     
    if ( abs(s2-s1)> e) then         
		s1 = s2
		s1n = s2n
                                 
    else                         
		s = s2                   
		exit
    end if      
end do
end subroutine trapezios
!-----------------------------------
subroutine cilindro(p,L,R,xc,d,xo,eps)
Real(prc), intent(in)   :: xc,p,L,R
Real(prc), intent(out)  :: d
Real(prc), intent(in):: xo, eps
a=xo-xc

call trapezios(funcao,-R,R,eps,d)

End subroutine   
!-----------------------------------
End module
program integral
use Grav_cilindro
implicit none
integer:: N, g,nc
real(prc):: Y0,fi
real(prc), allocatable:: x(:,:), y(:),ys(:), t(:)
real(prc), allocatable:: xtx(:,:), xty(:), pa(:)
integer:: i,j
real(prc):: eps, xc

open(22,file='dados_de_entrada.in',status='old',action='read')
open(33,file='grav_cilindros.dat',status='old',action='read') 

Read(22,*) N
Read(22,*) g
Read(22,*) nc
Read(22,*) eps
Read(22,*) R, L, p

allocate(x(n,nc),y(n),ys(n),xtx(nc,nc),xty(nc),pa(nc),t(n))

do i=1,n
   read(33,*) t(i),y(i)
end do
 close(33)

Do j=1,nc
   Read(22,*) xc
   Do i=1,n       
    !   a=t(i)-xc
      Call cilindro(p,L,R,xc,X(i,j),t(i),eps)       
   End Do    
End do

 call Monta_sistema_MQ(N, Nc, x, y, xtx, xty)
 call elimina(nc, xtx, xty)
 call retrosubs(nc, xtx, xty, pa)
 ys=matmul(x,pa)
 
write(*,*) pa(1), pa(2), pa(3) 
End program

