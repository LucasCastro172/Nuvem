program maior_valor
Implicit none
Real:: maior
Real:: x
Integer::l,i
Integer:: erro

Open(33,file='lista.dat',status='old',action='read')
read(33,*) maior
 
 l=1
 i=1
 
 Do
     read(33,*,iostat=erro) x    
     if(erro==0) then
        i = i + 1  
        if(x>maior) then
           maior=x
           L=i
        End if
     Else 
         exit     
     End if
     
 End do

write(*,*) 'O arquivo tem', i,'linhas'
Write(*,*) 'O maior valor é',maior,'na linha',L
End Program
 
