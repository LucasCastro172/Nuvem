%%%%%%Codigo Metodos Computacionais%%%%%%%%%%%%

%%%%%%%%Parametros%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numero de thread
th = 1:8
% parallel
par = [0.0428,0.0306,0.025,0.0279,0.027,0.0331,0.034,0.0303]
% serial
on = [0.0358,0.0576,0.0757,0.0769,0.0374,0.0332,0.0329,0.0312]

%%%%%%%%%%%Consideração de perfomance%%%%%%%%%%%
%%%%%%%%%Dividir thread 1 por todos os outros%%%
%%%%%S=t1/Tp
% parallel
spar = [1, 1.4,1.71,1.53,1.58,1.29,1.26,1.41]
% omp dentro
%son = [1, 1.99,1.60,2.07,1.69,1.73,1.78,1.99]

figure (1)
hold on
plot(th,par,'-r', 'LineWidth',2)
%plot(th,on,'--k','LineWidth',2)
%legend('Fora loop','No loop')
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Threads')
ylabel('Tempo (seg)' )
title('Threads x Tempo')

figure (2)
hold on
plot(th,spar,'-r', 'LineWidth',2)
%plot(th,son,'--k','LineWidth',2)
%legend('Fora loop','No loop')
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Threads')
ylabel('Speedup' )
title('Performance')