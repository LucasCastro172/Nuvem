!Programa para relacionar o consumo de gasolina por km andados
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021

program consumo_gasolina
Real:: l, km, con !parametros de kilometros, gasolina e consumo
real :: kmt,litro_total ! parametros para acumular a kilometragem e o gasto de gasolina
Real::soma,consumo !parametro para fazer a soma 
integer:: i !contador do loop
i=0.
soma=0.
km=0.
kmt=0.
litro_total=0.

Do
   i=i+1
   write(*,*)'Quantos litros?'
   Read(*,*)l
   If(l==-1)exit !número utilizado para encontrar o total 
     write(*,*)'Quantos quilômetros?'
     Read(*,*)km
     kmt=kmt+km
     con=km/l
     litro_total = litro_total + l
     Write(*,*)kml
     soma=soma+con
     kmt=kmt+km
     litro_total = litro_total + l
End do

consumo=kmt/litro_total
if(l == -1) then
write(*,*)'O consumo total foi de',consumo,'km/l'
end if


end program
