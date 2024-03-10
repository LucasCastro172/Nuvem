Program salario_vendedores 
character(9)::  nome, mes
Real:: preco
Real:: salario_B, salario_C, salario_A,preco_total_B,preco_total_C,preco_total_A
Integer:: i, n 

Open(33,file='AMANDA GUIMARAES PEREIRA - vendas_01.dat',status='old',action='read')
Open(44,file='AMANDA GUIMARAES PEREIRA - vendas_02.dat',status='old',action='read')
Open(55,file='AMANDA GUIMARAES PEREIRA - vendas_03.dat',status='old',action='read')

Read(33,*) mes
Read(33,*)
Read(33,*) N
Read(33,*)
Read(33,*)
Read(33,*)
Write(*,*) N
preco_total_B=0.
preco_total_C=0.
preco_total_A=0.

Do i=1,n
 Read(33,*)nome,preco
 call salario(n,nome,mes,preco,salario_B, salario_C, salario_A,preco_total_B,preco_total_C,preco_total_A)
End Do
write(*,*)'Salário de ',mes  
write(*,*) 'Alice', salario_A
write(*,*) 'Bruno', salario_B
write(*,*) 'Cristina', salario_C

End program
!=================
Subroutine salario(n,nome, mes,preco,salario_B, salario_C, salario_A,preco_total_B,preco_total_C,preco_total_A)
character(9),intent(in):: nome, mes
Real, intent(in):: preco
Real,intent(out):: salario_B, salario_C, salario_A
Real,intent(inout)::preco_total_B,preco_total_C,preco_total_A
Integer:: i
Integer,intent(in):: n 

 !Do i=1,N
   ! Read(33,*)nome,preco
   If(nome=='Bruno') then
       preco_total_B=preco_total_B+preco
       If(preco_total_B<60000) then
          salario_B=preco_total_B*5/100+1045
       Else if(preco_total_B>=60000.and.preco_total_B<70000) then
       salario_B=preco_total_B*7/100+1045
       Else if(preco_total_B>=70000) then
       salario_B=preco_total_B*7/100+1045+500
       End if
   Else If(nome=='Cristina') then

       preco_total_C=preco_total_C+preco       
       If(preco_total_C<60000) then
          salario_C=preco_total_C*5/100+1045
       Else if(preco_total_C>=60000.and.preco_total_C<70000) then
       salario_C=preco_total_C*7/100+1045
       Else if(preco_total_C>=70000) then
       salario_C=preco_total_C*7/100+1045+500
       End if
   Else If(nome=='Alice') then
       preco_total_A=preco_total_A+preco
       If(preco_total_A<60000) then
          salario_A=preco_total_A*5/100+1045
       Else if(preco_total_A>=60000.and.preco_total_A<70000) then
       salario_A=preco_total_A*7/100+1045
       Else if(preco_total_A>=70000) then
       salario_A=preco_total_A*7/100+1045+500
       End if
       End if
    ! End Do
End Subroutine
