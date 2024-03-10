Program converter_rad_graus
Real:: rad, graus, minu,seg
!Real, parameter:: pi=3.14159265359
Real, parameter:: pi=3.1415926
Write(*,*)'Valor em radianos?'
Read(*,*) rad

graus=rad*180/pi
minu=(graus-int(graus))*60
seg=(minu-int(minu))*60

write(*,*) rad,'rad=',int(graus),'graus,',int(minu),'min e', seg,'seg'	

End program  
