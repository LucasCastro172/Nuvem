clc
clear
close all;

dd=load('resultado_ex7.dat');
i=dd(:,1);
a_i=dd(:,2);
%a_i1=dd(:,7);
%dd1=load('dados_suavizado.dat')
plot(i,a_i,'b','DisplayName','Sol.Em','linewidth',1.5);
grid
