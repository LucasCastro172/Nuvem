module Sistemas
implicit none 
integer, parameter:: prc=kind(1.d0)

contains

subroutine retrosubs(N, U, B, x)
! A matriz 
implicit none 
integer, intent(in):: N
real(prc), intent(in):: U(N,N)
real(prc), intent(in):: b(N)
real(prc), intent(out):: x(N)
integer:: i,j
real(prc):: soma

x(N)=b(N)/U(N,N)

do i=N-1,1,-1
  Soma=0._prc
   do j=i+1,N
      Soma = Soma + U(i,j)*x(j)
   End Do
   x(i)=(b(i)- Soma)/U(i,i)
       
End do 
End Subroutine retrosubs
!!!!!!!!!!!!!!!!!!!!!!!
subroutine elimina(N,A,b)
! presupõem que a matriz é singular 
implicit none
integer,intent(in):: N
real(prc),intent(inout):: A(N,N),b(N)
real(prc):: m
integer:: i,j,k

Do k=1,n-1
   Do i=k+1,n
      m=A(i,k)/A(k,k)           
      Do j=k+1,n
         A(i,j)=A(i,j)-m*A(k,j)        
      End do
     ! A(i,k+1:N)=A(i,k+1:N)-m*A(k,k+1:N)
      b(i)=b(i)-m*b(k)
    End Do   
End Do

End subroutine elimina
!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine elimina_Pivoteamento(N,A,b,erro)
implicit none
integer,intent(in):: N
real(prc),intent(inout):: A(N,N),b(N)
real(prc):: m
integer:: i,j,k, I_max
integer,intent(out):: erro
real(prc):: maior,aux
real(prc):: eps=1.d-8
! Na forma de vetor:
real(prc):: Linha(n)

erro=0
Do k=1,n-1
   maior = abs(A(k,k))
   I_max=k
    do i = k+1, n
         if(abs(A(i,k))> Maior) then 
             maior = abs(A(i,k))
             I_max = i
         End if
    End do      
 If(maior > eps) then     
    linha(k:n) = A(k,k:n)
    A(k,k:n) = A(I_max,k:n)
    A(I_max,k:n)=linha(k:n)
    aux=b(k)
    b(k)=b(i_max)
    b(i_max)=aux

   Do i=k+1,n
      m = A(i,k) / A(k,k)           
      Do j=k+1,n
         A(i,j)=A(i,j)-m*A(k,j)        
      End do
     ! A(i,k+1:N)=A(i,k+1:N)-m*A(k,k+1:N)
      b(i)=b(i)-m*b(k)
    End Do   
  Else 
    erro=1
    write(*,*)'Matriz sigular ou quase.'
   exit
  End if  
End Do
End subroutine elimina_Pivoteamento 
Subroutine fatora_LU(N,A)
implicit none 
integer, intent(in):: N
real(prc):: A(N,N)
real(prc):: m
integer:: i,j,k

Do k=1,n-1
   Do i=k+1,n
      m=A(i,k)/A(k,k)           
      Do j=k+1,n
         A(i,j)=A(i,j)-m*A(k,j)        
      End do
      A(i,k)=m
    End Do   
End Do
End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine resolve_LU(N, LU, B)
! A matriz 
implicit none 
integer, intent(in):: N
real(prc), intent(in):: LU(N,N)
real(prc), intent(inout):: b(N)
integer:: i,j
real(prc):: soma

do i=2,N
   soma=0.d0
   do j=1,i-1
      soma=soma +LU(i,j)*b(j)
   end do 
   b(i) = b(i) - soma
end do       

b(N)=b(N)/LU(N,N)

do i=N-1,1,-1
  Soma=0._prc
   do j=i+1,N
      Soma = Soma + LU(i,j)*b(j)
   End Do
   b(i)=(b(i)- Soma)/LU(i,i)
       
End do 
End Subroutine resolve_LU
!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Monta_sistema_MQ_reta(N, x, y, xtx, xty) 
implicit none 
integer,intent(in):: N
real(prc), intent(in):: X(n,2), y(n)
real(prc), intent(out):: xtx(2,2), xty(2)
integer:: i,j,k

!xtx=matmul(transpose(x),x)
!xty=matmul(transpose(x),y)

do i=1,2
   do j=i,2
      xtx(i,j)=dot_product(x(:,i),x(:,j))
      xtx(j,i)=xtx(i,j)
   end do 
   xty(i)=dot_product(x(:,i),y)
end do 
      
End subroutine Monta_sistema_MQ_reta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Monta_sistema_MQ(N, Nc, x, y, xtx, xty) 
implicit none 
integer,intent(in):: N, Nc
real(prc), intent(in):: X(n,nc), y(n)
real(prc), intent(out):: xtx(nc,nc), xty(nc)
integer:: i,j,k

!xtx=matmul(transpose(x),x)
!xty=matmul(transpose(x),y)

do i=1,nc
   do j=i,nc
      xtx(i,j)=dot_product(x(:,i),x(:,j))
      xtx(j,i)=xtx(i,j)
   end do 
   xty(i)=dot_product(x(:,i),y)
end do 
      
End subroutine Monta_sistema_MQ
End module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!11!!!!!!!
!Program teste_U
!use Sistemas
!implicit none 
!integer:: N
!real(prc), allocatable:: U(:,:)
!real(prc), allocatable:: b(:), x(:),x0(:)
!integer:: i,j

!write(*,*)'N?'
!Read(*,*)N
!Allocate(U(N,N),x(N),b(N))

!U = 2*U + 1.d0
!  Do i=2,N
!    U(i,1:i-1)=0.d0
!  End do
  
!x0=[(i,i=1,N)]
!b=matmul(U,x0)
!  
!call retrosubs(N,U,b,x)
!
!do i=1,N
!  write(*,*)x0(i),x(i)
!End do
!End program 
!!!!!!!!!!!!!!!!!!!!!!!!!1
!Program teste_U
!use Sistemas
!implicit none 
!integer:: N, g,nc
!real(prc), allocatable:: x(:,:), y(:),ys(:),deltay(:)
!real(prc), allocatable:: xtx(:,:), xty(:), p(:)
!integer:: i

!n=100
!g=3
!nc=g+1 

!allocate(x(n,nc),y(n),ys(n),xtx(nc,nc),xty(nc),p(nc),deltay(n))

!x(:,1)=1.d0
!open(55, file='dados_para_mq.dat', status='old',action='read')
!do i=1,n
!   read(55,*) x(i,2), y(i)
!end do
! close(55)
 
! do i = 3, nc
!   x(:,i) = x(:,2)**(i-1)
! end do        
! call Monta_sistema_MQ(N, Nc, x, y, xtx, xty)
 !call fatora_LU(n,xtx)
 !call resolve_LU(N,xtx,xty)
! call elimina(nc, xtx, xty)
! call retrosubs(nc, xtx, xty, p)
! write(*,*) 'A:', p(1)
! write(*,*) 'B:', p(2)
! write(*,*) 'C:', p(3)
! write(*,*) 'D:', p(4)
! ys=matmul(x,p)
!open(22,file='P2.out',status='replace')
!do i=1,N
!  deltay(i) = ys(i) - y(i)
!  write(*,*) x(i,2), y(i), ys(i)
!End do
!do i=1,N
!  write(22,*) x(i,2), y(i), ys(i),deltay(i)
!End do

!End program 
