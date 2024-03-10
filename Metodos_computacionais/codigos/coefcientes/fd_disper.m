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
MM=4*M;
A=zeros(MM,N);
b=zeros(MM,1);
L=zeros(N+2,N+2);
bp= zeros(N+2,1);
phi=[0 pi/12 pi/6 pi/4.0];
%
% matching spectral properties
%
k=0.5*[0:M-1]'*pi/M;
%b=-k.^2;
M0=-M+1;
M1=0;
for ip=1:4
M0=M0+M;
M1=M1+M;
cphi=cos(phi(ip));
sphi=sin(phi(ip));
Wa = cphi*cphi*sqrt(1+(k/max(k)).^2);
%Wa = cphi^2;
b(M0:M1)=-2.0*diag(Wa)*(sin(k*mu/2)/mu).^2;
A(M0:M1,1) = Wa;
A(M0:M1,2:N) = diag(Wa)*(cos(cphi*k(1:M)*[1:N-1])+cos(sphi*k(1:M)*[1:N-1]));
end
%
L( 1 ,1:N) = A(1:MM,1)'*A(1:MM,1:N);   L(1,N+2) =-1.0; 
L(2:N,1:N) = A(1:MM,2:N)'*A(1:MM,1:N); L(2:N,N+1)=-[1:N-1]'.^2; L(2:N,N+2) =-2.0; 
L(N+1,2:N)=[1:N-1].^2; 
L(N+2,1)   = 1.0; L(N+2,2:N)=2.0; 
%
bp(1:N,1)=A'*b;
bp(N+1,1)=1.0;
%
clsc = L \ bp;
%
% Matching phase and group
%
MM=8*M;
A=zeros(MM,N);
b=zeros(MM,1);
L=zeros(N+2,N+2);
bp= zeros(N+2,1);
ein=ones(M,1);
%
k=0.5*[0:M-1]'*pi/M;
%b=-k.^2;
M0=-M+1;
M1=0;

for ip=1:4
M0=M0+M;
M1=M1+M;
cphi=cos(phi(ip));
sphi=sin(phi(ip));
Wa = cphi*cphi*sqrt(1+(k/max(k)).^2);
b(M0:M1)     =-2.0*diag(Wa)*(sin(k*mu/2)/mu).^2;
A(M0:M1,1)   = Wa;
A(M0:M1,2:N) = diag(Wa)*(cos(cphi*k(1:M)*[1:N-1])+cos(sphi*k(1:M)*[1:N-1]));
%
M0=M0+M;
M1=M1+M;
Wa = cphi;
b(M0:M1)     = diag(Wa)*(sin(k*mu)/mu);
A(M0:M1,1)   = 0.0;
A(M0:M1,2:N) = diag(Wa)*(cphi*sin(cphi*k(1:M)*[1:N-1])+sphi*sin(sphi*k(1:M)*[1:N-1])).*(ein*[1:N-1]);
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
% Matching Polynomial derivatives 
%
A=zeros(N,N);
b=zeros(N,1);
b(2)=1.0;
A(1:1) = 0.5;
A(1,2:N) = 1.0;
for m=2:N
    A(m,2:N) = ([1:N-1]).^(2*(m-1)); 
end
c = A \ b;
%
% Dispersion analysis
%
Nd=360;
wn=[1:Nd]'*pi/(Nd+1);
kn=[1:Nd]'*pi/(Nd+1);
cfd=0.0*kn;vg=0.0*kn;
ein=ones(Nd,1);

d = c;
for ip=1:4
cphi=cos(phi(ip));
sphi=sin(phi(ip));
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(d(1)*ein + ...
                    cos(kn*cphi*[1:N-1])*d(2:N) + ...
                    cos(kn*sphi*[1:N-1])*d(2:N) ))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*cphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*sphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 ) ./ sin(mu*cfd.*kn);

figure(2),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
set(h,'LineWidth',2);
if ip == 1
  set(gca,'FontWeight','bold','FontSize',14);
  xlabel('dx/\lambda');
  legend('c fase','V grupo')
  hold on
end 
end
hold off;
title('Derivada de polinomios')
%
%
%
d = clsc(1:N);

for ip=1:4
cphi=cos(phi(ip));
sphi=sin(phi(ip));
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(d(1)*ein + ...
                    cos(kn*cphi*[1:N-1])*d(2:N) + ...
                    cos(kn*sphi*[1:N-1])*d(2:N) ))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*cphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*sphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 ) ./ sin(mu*cfd.*kn);
figure(3),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
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
%
%
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
scale=sqrt(0.5*abs(d(1))/dot(d1,d1));
d1=scale*d1;


for ip=1:4
cphi=cos(phi(ip));
sphi=sin(phi(ip));
%
% phase velocity
%
cfd = (2.0/mu)*asin(sqrt(-0.5*mu*mu*(d(1)*ein + ...
                    cos(kn*cphi*[1:N-1])*d(2:N) + ...
                    cos(kn*sphi*[1:N-1])*d(2:N) ))) ./ kn;
%
% group velocity
%
vg = mu*sqrt((sin(kn*cphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 + ...
             (sin(kn*sphi*[1:N-1])*(d(2:N).*[1:N-1]')).^2 ) ./ sin(mu*cfd.*kn);
figure(4),h=plot(kn/pi/2,cfd,'r',kn/pi/2,vg,'.b',kn/pi/2,ein,'k');axis([0 0.5 0.90 1.10]);
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
d1clscv=sqrt(0.5*abs(d(1))/dot(d1,d1))*d1;


coef=[1 2.0*cos([1:N-1]*pi*sqrt(3.0)/3.0)];
disper2d= sqrt(-2.0/dot(coef,clscv(1:N)))

mu/disper2d

