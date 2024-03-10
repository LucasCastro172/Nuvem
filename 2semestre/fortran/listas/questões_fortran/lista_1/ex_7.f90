program consumo_gasolina
Real:: l,km,kml,kmt,litro_total
Real::soma,consumo
character(1) SouN
integer:: i
i=0.
soma=0.
km=0.
kmt=0.
litro_total=0.

Do
   i=i+1
   write(*,*)'Quantos litros?'
   Read(*,*)l
   If(l==-1)exit
     write(*,*)'Quantos quilômetros?'
     Read(*,*)km
     kmt=kmt+km
     kml=km/l
     litro_total = litro_total + l
     Write(*,*)kml
!write(*,'(A,x,F3.5,x,A)')'Consumo:',kml,'km/l'
     soma=soma+kml
     kmt=kmt+km
     litro_total = litro_total + l
     Write(*,*)soma
!write(*,*)'Deseja continuar?(S ou N)'
!Read(*,*)SouN
!If(SouN=='N'.or.SouN=='n') exit
End do

consumo=kmt/litro_total

write(*,'(A,F15.5,x,A)')'Consumo:',consumo,'km/l'

End program
