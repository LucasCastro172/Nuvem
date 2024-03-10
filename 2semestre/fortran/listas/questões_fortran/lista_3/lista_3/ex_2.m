clc
clear
close all;
Ax = load('P1.out');
i = Ax(:,1) ;
a_i = Ax(:,2) ;
plot(i,a_i,'b','DisplayName','diferenças finitas','linewidth',1.7);
hold on
i = Ax(:,1) ;
a_i1 = Ax(:,3) ;
plot(i,a_i1,'r','DisplayName','diferenças finitas','linewidth',1.7);
hold on
Ax1 = load('P2.out');
i = Ax1(:,1) ;
a_i3 = Ax1(:,3) ;
plot(i,a_i3,'k','DisplayName','diferenças finitas','linewidth',1.7);
grid

%xticks(-50:10:50)
%yticks(0:10:60)
%axis([-50, 50, 0, 60]);
%title('Gravimetria do cilindro','Interpreter','latex','FontSize', 22);
%xlabel('a(m)','FontSize', 16);
%ylabel('g/Gp(m.s)','FontSize', 16)
