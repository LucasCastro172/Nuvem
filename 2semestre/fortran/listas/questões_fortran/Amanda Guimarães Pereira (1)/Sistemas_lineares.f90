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
