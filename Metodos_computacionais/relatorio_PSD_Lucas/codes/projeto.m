%%%%%%%%Projeto para o curso de "Processamento de sinais digitais"%%%%%%%%%%%%%%%
%%%%%%%%%%%%Multi convolucoes: funcao retangular%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###########################
%%%%%%%%%Parametros utilizados
Fs = 1000;              %parametro relacionado ao intervalo de amostragem do tempo
dt = 1/Fs;              %equacao para ointervalo de amostragem do tempo
%delf  =  1/(n*dt)      %equacao para o intervalo de amostragem da frequencia

################################################################################
%%%%%%%%%%%%%%%%%%%%%%%%%%%Funcao caixa retangular%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rect=@(x,a) ones(1,numel(x)).*(abs(x)<a/2); % a define a largura do pulso
x=-5:dt:5;          %Vetor para o numero de pontos 
y=rect(x,2);        %Formato da funcao retangular, escolhendo a=2
ynorm=y/(max(y));   %Normalizacao da funcao retangular

 ###################Plotagem funcao retangular##################################
figure (99) 
plot(x,ynorm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Funcao retangular')
%%%%%%%%%  Espectro de amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(ynorm);     %Aplicacao da transformada rapida de Fourier na funcao retangular
L = 10000;          %Tamanho do eixo x, numero de pontos
P2 = abs(Y/L);      %Formula para a amplitude 

%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1norm = P1/max(P1);          %Normalizacao do resultado
f = Fs*(0:(L/2))/L;

%%%%%%%%%%%%%Plotagem do espectro de amplitude da funcao retangular
figure(999)
plot(f,P1norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 5])
xlabel('Numero de pontos')
ylabel('Amplitude')
#######################################################
#######################################################
%Primeira convolucao
c1=conv(y,y)*dt;     %Aplicacao da auto convolucao 
c1norm=c1/(max(c1)); %Normalizacao do convolucao 
x1=-10:dt:10;        %Vetor com o numero de pontos
%%%%%%%%%Plotagem da convolucao 1
figure (1)
plot(x1,c1norm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Primeira convolucao')

%%%%%%%%%%%%%%%Espectro de amplitude
c1ft = fft(c1norm);  %Aplicacao da transformada rapida de Fourier
L1 = 20000;          %Tamanho do eixo x, numero de pontos
P2c1 = abs(c1ft/L1); %calculo para a amplitude
%df1=0.005           % Valor do intervalo de amostragem na frequencia

%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1c1 = P2c1(1:L1/2+1);
P1c1(2:end-1) = 2*P1c1(2:end-1);
P1c1norm = P1c1/max(P1c1);%Normalizacao do resultado
f1 = Fs*(0:(L1/2))/L1;

%%%%%%%%%%%%%%Plotagem do espectro de amplitude
figure(100)
plot(f1,P1c1norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 2])
xlabel('Numero de pontos')
ylabel('Amplitude')
#######################################################
#######################################################
%Segunda convolucao
c2=conv(c1,c1)*dt;   %Aplicacao da auto convolucao do resultado anterior
c2norm=c2/(max(c2)); %Normalizacao do resultado
x2=-20:dt:20;        %Vetor do eixo x, numero de pontos

%%%%%%%%%%%Plotagem da convolucao 2
figure (2)
plot(x2,c2norm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Segunda convolucao')
%%%%%Espectro de amplitude

c2ft = fft(c2norm); %Aplicacao da transformada rapida de Fourier
L2 = 40000;         %Tamanho do eixo x, numero de pontos
P2c2 = abs(c2ft/L2);%calculo para a amplitude
%df2=0.025          % Valor do intervalo de amostragem na frequencia
%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1c2 = P2c2(1:L2/2+1);
P1c2(2:end-1) = 2*P1c2(2:end-1);
P1c2norm = P1c2/max(P1c2);       %Normalizacao do resultado
f2 = Fs*(0:(L2/2))/L2;

%%%%%%%%%%%%%%Plotagem do espectro de amplitude
figure(200)
plot(f2,P1c2norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 2])
xlabel('Numero de pontos')
ylabel('Amplitude')

##################################
%Terceira convolucao
%%%%numero de pontos
c3=conv(c2,c2)*dt;      %Aplicacao da auto convolucao do resultado anterior
c3norm=c3/(max(c3));    %Normalizacao do resultado
x3=-40:dt:40;           %Vetor do eixo x, numero de pontos

%%%%%%%%%%%Plotagem da convolucao 3
figure (3)
plot(x3,c3norm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Terceira convolucao')
%Espectro de amplitude
c3ft = fft(c3norm);    %Aplicacao da transformada rapida de Fourier
L3 = 80000;            %Tamanho do eixo x, numero de pontos
P2c3 = abs(c3ft/L3);   %calculo para a amplitude
%df3=0.0125            %Valor do intervalo de amostragem na frequencia

%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1c3 = P2c3(1:L3/2+1);
P1c3(2:end-1) = 2*P1c3(2:end-1);
P1c3norm = P1c3/max(P1c3);%Normalizacao do resultado
f3 = Fs*(0:(L3/2))/L3;

%%%%%%%%%%%%%%Plotagem do espectro de amplitude
figure(300)
plot(f3,P1c3norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 2])
xlabel('Numero de pontos')
ylabel('Amplitude')

######################################
%Quarta convolucao 
c4=conv(c3,c3)*dt;     %Aplicacao da auto convolucao do resultado anterior
c4norm=c4/(max(c4));   %Normalizacao do resultado
x4=-80:dt:80.;         %Vetor do eixo x, numero de pontos

%%%%%%%%%%%Plotagem da convolucao 2
figure (4)
plot(x4,c4norm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Quarta convolucao')

%%%%%%%%Amplitude espectral
c4ft = fft(c4norm);    %Aplicacao da transformada rapida de Fourier
L4 = 160000;           %Tamanho do eixo x, numero de pontos
P2c4 = abs(c4ft/L4);   %calculo para a amplitude
%df4=0.00625           %Valor do intervalo de amostragem na frequencia
%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1c4 = P2c4(1:L4/2+1);
P1c4(2:end-1) = 2*P1c4(2:end-1);
P1c4norm = P1c4/max(P1c4);       %Normalizacao do resultado
f4 = Fs*(0:(L4/2))/L4;

%%%%%%%%%%%%%%Plotagem do espectro de amplitude
figure(400)
plot(f4,P1c4norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 2])
xlabel('Numero de pontos')
ylabel('Amplitude')

################################################################
%quinta convolucao
c5=conv(c4,c4)*dt;   %Aplicacao da auto convolucao do resultado anterior
c5max= max(c5);
c5norm= c5/c5max;    %Normalizacao do resultado
x5=-160:0.001:160.;  %Vetor do eixo x, numero de pontos

%%%%%%%%%%%Plotagem da convolucao 2
figure (5)
plot(x5,c5norm,'LineWidth',2)
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'Xlim',[-6 6])
xlabel('numero de pontos')
ylabel('Amplitude' )
title('Quinta convolucao')

%Amplitude espectral
c5ft = fft(c5norm);   %Aplicacao da transformada rapida de Fourier
L5 = 320000;          %Tamanho do eixo x, numero de pontos
P2c5 = abs(c5ft/L5);  %calculo para a amplitude
%df5=0.003125         %Valor do intervalo de amostragem na frequencia


%%%%%%%%%%%%Calculo feito para considerar apenas uma banda da amplitude
P1c5 = P2c5(1:L5/2+1);
P1c5(2:end-1) = 2*P1c5(2:end-1);
P1c5norm = P1c5/max(P1c5);        %Normalizacao do resultado
f5 = Fs*(0:(L5/2))/L5;

%%%%%%%%%%%%%%Plotagem do espectro de amplitude
figure(500)
plot(f5,P1c5norm,'LineWidth',2) 
title('Espectro de amplitude ')
set(gca,'Xlim',[0 2])
xlabel('numero de pontos')
ylabel('Amplitude')

###################################################################################
%%%%%%%%%%%%graficos comparativos%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Comparacao das convolucoes e a funcao normalizados%%%%%%
figure(101)
hold on
plot(x,ynorm,'-k','LineWidth',2)
plot(x1,c1norm,'-y','LineWidth',2)
plot(x2,c2norm,':c', 'LineWidth',2)
plot(x3,c3norm,'-.m', 'LineWidth',2)
plot(x4,c4norm,'-r','LineWidth',2)
plot(x5,c5norm,'--b','LineWidth',2)
set(gca,'Xlim',[-6 6])
legend('Funcao','Primeira','Segunda','Terceira','Quarta','Quinta' )
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Numero de pontos')
ylabel('Amplitude' )
title('Convolucoes')

%%%%%%%%%%%%%Comparacao Amplitude Espectral das convolucoes e da funcao normalizada%%

figure(102)
hold on

plot(f,P1norm,'-k','LineWidth',2)
plot(f1,P1c1norm,'-y','LineWidth',2)
plot(f2,P1c2norm,':c', 'LineWidth',2)
plot(f3,P1c3norm,'-.m', 'LineWidth',2)
plot(f4,P1c4norm,'-r','LineWidth',2)
plot(f5,P1c5norm,'--b','LineWidth',2)
set(gca,'Xlim',[0 2])
legend('Funcao','Primeira','Segunda','Terceira','Quarta','Quinta' )

xlabel('Numero de  pontos')
ylabel('Amplitude' )
title('Amplitude espectral')

