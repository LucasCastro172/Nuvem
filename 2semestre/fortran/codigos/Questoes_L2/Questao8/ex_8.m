clc
clear
close all;

dd=load('resultado_ex8.dat');
i=dd(:,1);
a_i=dd(:,2);
%a_i1=dd(:,7);
%dd1=load('dados_suavizado.dat')
i_1=dd(:,3);
a_i1=dd(:,4);
f1 = figure(1);
plot(i,a_i,'b','DisplayName','Sol.Em','linewidth',1.5);
hold on
plot(i_1,a_i1,'r','DisplayName','Sol.Analitica','linewidth',1.5);
grid
