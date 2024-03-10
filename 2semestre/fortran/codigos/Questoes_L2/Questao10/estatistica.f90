module estatistica 
implicit none 
Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
Real function media(N,v)
implicit none 
integer, intent(in):: N
Real, intent(in):: V(n)

media=sum(v)/N

End function media
!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine md_dp(N,v,md,dp)
implicit none
integer, intent(in):: N
Real, intent(in):: V(n)
Real, intent(out):: md,dp
Real:: soma
Integer:: i

md=media(n,v)

dp=sqrt(sum((v-md)**2)/n)
!soma=0.
!Do i=1,n
!   soma=soma+(v(i)-md)**2
!End Do   
!   dp=sqrt(soma/N)
   
End subroutine md_dp
!!!!!!!!!!!!!!!!!!!!!!!
Subroutine ordena(N,v)
implicit none 
integer, intent(in):: N
real:: V(N)
integer:: i,j,k, troca
real:: aux ! variavel temporaria para guardar o valor durante as trocas

k=0
 Do i=1,N-1
    
    troca=0
    Do j=1,N-i
       if(v(j)>v(j+1)) then
          aux=v(j+1)
          v(j+1)=v(j)
          v(j)=aux
          troca=1
       End if       
    End do
    k=k+1
    if(troca==0) exit
 
 End Do   
 
 End subroutine ordena 
!!!!!!!!!!!!!!!!!!!!!!!
Subroutine mediana(N,v,mn)
implicit none 
integer, intent(in):: N
Real:: V(n)
real:: mn

 If(mod(N,2)==0) then
   mn= (v(N/2)+ v(N/2+1))/2
 Else
   mn= v((N-1)/2+1) 
 End if
 
End subroutine mediana 
End module estatistica
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program teste_vetores
!use estatistica
!implicit none

!Real, allocatable:: x(:)
!integer:: N
!Real:: mn
!write(*,*)'Quantos números?'
!read(*,*)N
!Allocate(x(n))

!call random_number(x)
!x=10*(x-0.5)

!write(*,*)x

!call ordena(N,x)

!call mediana(n,x,mn)

!Open(44,file='vetor_ordenado.dat',status='replace')
! write(44,'(f15.4)') mn
! write(44,'(f15.4)') x
! close(44)

!End program 

!Program teste_vetores
!use estatistica
!implicit none

!Real, allocatable:: x(:)
!integer:: N
!integer:: i
!Real:: m, dp
!Open(44,file='vetor.dat',status='old', action='read')
!Read(44,*) N
!Allocate(x(n))
!Do i=1,n
!   Read(44,*)x(i)
!End Do   

! close(44)
!!!!!!!!!!!!
!call md_dp(N,x,m,dp)
! m=media(n,x)
! write(*,*)'media:', m
! write(*,*)'Desvio padrão:',dp
!End program 
