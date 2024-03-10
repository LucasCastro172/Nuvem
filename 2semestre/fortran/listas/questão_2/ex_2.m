clc
clear
close all;
Ax = load('resultado_ex2.dat');
i = Ax(:,1) ;
a_i = Ax(:,2) ;
plot(i,a_i,'b','DisplayName','diferenças finitas','linewidth',1.7);
grid
xticks(-50:10:50)
yticks(0:10:60)
axis([-50, 50, 0, 60]);
title('Gravimetria do cilindro','Interpreter','latex','FontSize', 22);
xlabel('a(m)','FontSize', 16);
ylabel('g/Gp(m.s)','FontSize', 16)