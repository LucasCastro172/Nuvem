%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CALCULO DO OPERADOR FD:
%
% Dispersion relation for FD
% approximation of acoustic wave equation
%
% MINI-CURSO:
% Semana de Inverno UNICAMP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fitting the spectrum of the operator
%
global d;
N=7;          % number of independent coefficients
M=256;        % number of samples from spectrum
mu=0.20;
alpha=0.575;

A=zeros(M,N);
b=zeros(M,1);
L=zeros(N+2,N+2);
bp= zeros(N+2,1);
%
% matching spectral properties
%
phi=0.0;
cphi=cos(phi);
sphi=sin(phi);
k=alpha*[0:M-1]'*pi/M;
Wk = sqrt(1.0 + sqrt((k/max(k)).^2));
b=-k.^2;
A(1:M,1) = 1.0;
A(1:M,2:N) = 2.0*cos(k(1:M)*[1:N-1]);
%
L( 1 ,1:N) = A(1:M,1)'*A(1:M,1:N);   L(1,N+2) =-1.0; 
L(2:N,1:N) = A(1:M,2:N)'*A(1:M,1:N); L(2:N,N+1)=-[1:N-1]'.^2; L(2:N,N+2) =-2.0; 
L(N+1,2:N)=[1:N-1].^2; 
L(N+2,1)   = 1.0; L(N+2,2:N)=2.0; 
%
bp(1:N,1)=A'*b;
bp(N+1,1)=1.0;
%
clsk = L \ bp;
%
figure(1),h=plot(k/pi/2,-b,k/pi/2,-A*clsk(1:N),'.r');
set(h,'LineWidth',2);
set(gca,'FontWeight','bold','FontSize',14);
xlabel('Numero de onda');
ylabel('Espectro do operador')
legend('Exato','FD')
%
% Matching fase velocity
%
tht=[0 pi/12 pi/6.0 pi/4.0 pi/2.0];
phi=[0 pi/12 pi/6.0 pi/4.0 pi/2.0];
Ndir=length(tht)*length(phi);
nu1=zeros(1,Ndir);
nu2=zeros(1,Ndir);
nu3=zeros(1,Ndir);
ip=0
for itht=1:length(tht)
for iphi=1:length(phi)
    ip=ip+1;
    nu1(ip)=sin(tht(itht))*cos(phi(iphi));
    nu2(ip)=sin(tht(itht))*sin(phi(iphi));
    nu3(ip)=cos(tht(itht));
    
end
end
MM=Ndir*M;
A=zeros(MM,N);
b=zeros(MM,1);
L=zeros(N+2,N+2);
bp= zeros(N+2,1);
%
% matching spectral properties
%
k=alpha*[0:M-1]'*pi/M;
%b=-k.^2;
M0=-M+1;
M1=0;
for ip=1:Ndir
M0=M0+M;
M1=M1+M;
Wa = Wk;
%Wa = sqrt(1+(k/max(k)).^2);
%Wa = cphi^2;
b(M0:M1)=-2.0*diag(Wa)*(sin(k*mu/2)/mu).^2;
A(M0:M1,1) = 1.5*Wa;
A(M0:M1,2:N) = ...
    diag(Wa)*(cos(nu1(ip)*k(1:M)*[1:N-1])+cos(nu2(ip)*k(1:M)*[1:N-1])+cos(nu3(ip)*k(1:M)*[1:N-1]));
end
%
L( 1 ,1:N) = A(1:MM,1)'*A(1:MM,1:N);   L(1,N+2) =-1.0; 
L(2:N,1:N) = A(1:MM,2:N)'*A(1:MM,1:N); L(2:N,N+1)=-[1:N-1]'.^2; L(2:N,N+2) =-2.0; 
L(N+1,2:N) = [1:N-1].^2; 
L(N+2,1)   = 1.0; L(N+2,2:N)=2.0; 
%
bp(1:N,1)=A'*b;
bp(N+1,1)=1.0;
%
clsc = L \ bp;
%
% Matching phase and group
%
MM=2*Ndir*M;
A=zeros(MM,N);
b=zeros(MM,1);
L=zeros(N+2,N+2);
bp= zeros(N+2,1);
ein=ones(M,1);
%
k=alpha*[0:M-1]'*pi/M;
%b=-k.^2;
M0=-M+1;
M1=0;

for ip=1:Ndir
M0=M0+M;
M1=M1+M;
%Wa = Wk*sqrt(nu3(ip));
Wa = Wk;
b(M0:M1)     =-2.0*diag(Wa)*(sin(k*mu/2)/mu).^2;
A(M0:M1,1)   = 1.5*Wa;
A(M0:M1,2:N) = diag(Wa)*(cos(nu1(ip)*k(1:M)*[1:N-1]) ...
+cos(nu2(ip)*k(1:M)*[1:N-1]) ...
+cos(nu3(ip)*k(1:M)*[1:N-1]));
%
M0=M0+M;
M1=M1+M;
%Wa = 4.0*Wk*sqrt(nu3(ip));
Wa = 16.0*Wk;
b(M0:M1)     = diag(Wa)*(sin(k*mu)/mu);
A(M0:M1,1)   = 0.0;
A(M0:M1,2:N) = diag(Wa)*(nu1(ip)*sin(nu1(ip)*k(1:M)*[1:N-1])+...
                         nu2(ip)*sin(nu2(ip)*k(1:M)*[1:N-1])+...
                         nu3(ip)*sin(nu3(ip)*k(1:M)*[1:N-1])).*(ein*[1:N-1]);
end
%
L( 1 ,1:N) = A(1:MM,1)'*A(1:MM,1:N);   L(1,N+2) =-1.0; 
L(2:N,1:N) = A(1:MM,2:N)'*A(1:MM,1:N); L(2:N,N+1)=-[1:N-1]'.^2; L(2:N,N+2) =-2.0; 
L(N+1,2:N) = [1:N-1].^2; 
L(N+2,1)   = 1.0; L(N+2,2:N)=2.0; 
%
bp(1:N,1)=A'*b;
bp(N+1,1)=1.0;
%
clscv = L \ bp;
%
d = clsc(1:N);

%
% display results
%
Nd=360;
wn=[1:Nd]'*pi/(Nd+1);
kn=[1:Nd]'*pi/(Nd+1);
cfd=0.0*kn;vg=0.0*kn;
ein=ones(Nd,1);


for ip=1:Ndir
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(1.5*d(1)*ein + ...
                    cos(kn*nu1(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu2(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu3(ip)*[1:N-1])*d(2:N)))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*nu1(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu2(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu3(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2) ./ sin(mu*cfd.*kn);
figure(2),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
set(h,'LineWidth',2);

if ip == 1
  set(gca,'FontWeight','bold','FontSize',14);
  xlabel('dx/\lambda');
  legend('c fase','V grupo')
  hold on
end 
end
title('Ajuste Velocidade Fase')
hold off;
%
d = clscv(1:N);
cfd=0.0*kn;vg=0.0*kn;
for ip=1:Ndir
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(1.5*d(1)*ein + ...
                    cos(kn*nu1(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu2(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu3(ip)*[1:N-1])*d(2:N)))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*nu1(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu2(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu3(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2) ./ sin(mu*cfd.*kn);
figure(3),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg./cfd,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
set(h,'LineWidth',2);

if ip == 1
  set(gca,'FontWeight','bold','FontSize',14);
  xlabel('dx/\lambda');
  legend('c fase','V grupo')
  hold on
end 
end
title('Ajuste Velocidade Fase e Grupo')
hold off;

mu=0.30;

for ip=1:Ndir
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(1.5*d(1)*ein + ...
                    cos(kn*nu1(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu2(ip)*[1:N-1])*d(2:N) + ...
                    cos(kn*nu3(ip)*[1:N-1])*d(2:N)))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*nu1(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu2(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*nu3(ip)*[1:N-1])*(d(2:N).*[1:N-1]')).^2) ./ sin(mu*cfd.*kn);
figure(4),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg./cfd,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
set(h,'LineWidth',2);

if ip == 1
  set(gca,'FontWeight','bold','FontSize',14);
  xlabel('dx/\lambda');
  legend('c fase','V grupo')
  hold on
end 
end
title('Ajuste Velocidade Fase e Grupo D1')
hold off;

coef=[1.5 3.0*cos([1:N-1]*pi*sqrt(3.0)/3.0)];
disper3d= sqrt(-2.0/dot(coef,clscv(1:N)))

mu/disper3d

d = clsc(1:N);
d1=0*d;
d1(N)=sqrt(abs(d(N)));
d1(N-1)=-0.5*d(N-1)/d1(N);
for j=2:N-1
    aux=0.0;
    for iconv=1:j-1
        aux=aux+d1(N-iconv)*d1(N+iconv-j);
    end
    d1(N-j)=-0.5*(d(N-j)+aux)/d1(N);
end
d1clsc=sqrt(0.5*abs(d(1))/dot(d1,d1))*d1;

d = clscv(1:N);
d1=0*d;
d1(N)=sqrt(abs(d(N)));
d1(N-1)=-0.5*d(N-1)/d1(N);
for j=2:N-1
    aux=0.0;
    for iconv=1:j-1
        aux=aux+d1(N-iconv)*d1(N+iconv-j);
    end
    d1(N-j)=-0.5*(d(N-j)+aux)/d1(N);
end
d1clscv=sqrt(0.5*abs(d(1))/dot(d1,d1))*d1


