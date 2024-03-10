module matrizes
implicit none 

Contains 

subroutine produto_matriz_vetor(m,n,a,x,y)
implicit none 
integer, intent(in):: m,n
real, intent(in):: A(m,n), x(n)
real, intent(out):: y(m)
integer:: i, j

 Do i=1,m
   y(i)=dot_product(A(i,:),x)
 End Do 
 
 !Do i=1,m
 !   y(i)=0.
 !   Do j=1,N
 !      y(i)=y(i)+ A(i,j)*x(j)
 !   End do 
 !End do      
 
 End subroutine 
 End module matrizes 
 !!!!!!!!!!!!!!!!!!!!!
 program testaproduto
 use matrizes  
 implicit none 
 integer:: m,n
 real, allocatable:: A(:,:),x(:),y(:)
 integer:: i,j 
  call random_seed()
 write(*,*)'Números de linhas e colunas:'
 Read(*,*)m,n
 allocate(A(m,n),x(n),y(m))
 
 call random_number(A)
 A=int(20*(A-0.5))
  
  call random_number(x)
  x=int(10*(x-0.5))
  
  call produto_matriz_vetor(M,N,A,x,y)
  
  write(*,*)'Matriz A:'
    Do i=1,M
       write(*,*)(A(i,j),j=1,n)
    End do 
    
    write(*,*)'vetor x:'
    Do i=1,n
      write(*,*)x(i)
    End do
    
    write(*,*)'vetor y=Ax:'
    Do i=1,m
      write(*,*)y(i)
    End do 
End program      
          
