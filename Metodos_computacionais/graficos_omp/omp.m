%%%%%%%Codigo Metodos Computacionais%%%%%%%%%%%%

%%%%%%%%Parametros%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numero de thread
th = 1:8
% omp fora
off = [26.19,13.21,16.38,12.55,13.18,13.00,11.23,9.90]
% omp dentro
on = [26.00,13.09,16.29,12.53,15.41,15.02,14.58,13.07]

%%%%%%%%%%%Consideração de perfomance%%%%%%%%%%%
%%%%%%%%%Dividir thread 1 por todos os outros%%%
%%%%%S=t1/Tp
% omp fora
soff = [1, 1.98,1.60,2.09,1.99,2.01,2.33,2.64]
% omp dentro
son = [1, 1.99,1.60,2.07,1.69,1.73,1.78,1.99]

figure (1)
hold on
plot(th,off,'-r', 'LineWidth',2)
plot(th,on,'--k','LineWidth',2)
legend('Fora loop','No loop')
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Threads')
ylabel('Tempo (seg)' )
title('Comparação Fora x Dentro loop OMP')

figure (2)
hold on
plot(th,soff,'-r', 'LineWidth',2)
plot(th,son,'--k','LineWidth',2)
legend('Fora loop','No loop')
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Threads')
ylabel('Speedup' )
title('Performance')