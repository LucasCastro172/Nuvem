%
% Pos-graduacao em Geofisica - PPGGF (UFPA)
% Metodos Computacionais Aplicados a Sismica
% Prof. Jesse Costa
% 
% Deivid dos Santos Nascimento
% 30/09/2020
%

%
% Calculo dos coeficientes usados na aproximacao 
% de derivadas por Diferencas Finitas
%

% ordem da aproximacao
od = 8; %2; %4; %6;

n = od/2;

fprintf('Calculo de derivadas por diferencas finitas\n');
fprintf(' Coeficientes para aproximacao de %da ordem\n',od);
fprintf('===========================================\n');


%% coeficientes para primeira derivada (grid não intercalado)

%
% Coeficientes sao obtidos resolvendo-se o seguinte sistema:
% 
% kδ(k,1) = 2 Σ dj (2j + 1)^k   k = 1,3,...,2N-1
%    (somatorio de 0 a N-1)
%

%                                           %
% Monta e resolve o sistema 'Ax = b' -------%
%                                           %

% vetor k (1,3,...,2N-1)
k = 1:2:2*n-1;
nk = length(k);

% termo independente
b = zeros(nk,1);
b(1) = 1;

% matriz 'Akj'
A = zeros(nk,n);
for ik = 1:nk
    for ij = 0:n-1
        A(ik,ij+1) = ( ij + 1 )^k(ik);
    end
end
A = 2*A;

% resolve o sistema
d = linsolve(A,b);

%-------------------------------------------%


% Display do vetor de coeficientes 'd'
fprintf('     1a derivada\n');
fprintf('(grid não intercalado)\n');
fprintf('d=\n');% [');
for ik = 1:nk
    fprintf(' %14E\n',d(ik));
end
fprintf('--------------------\n');


%return


%% coeficientes para primeira derivada (grid intercalado)

%
% Coeficientes sao obtidos resolvendo-se o seguinte sistema:
% 
% kδ(k,1) = 2 Σ dj [(2j + 1)/2]^k   k = 1,3,...,2N-1
%    (somatorio de 0 a N-1)
%

%                                           %
% Monta e resolve o sistema 'Ax = b' -------%
%                                           %

% vetor k (1,3,...,2N-1)
k = 1:2:2*n-1;
nk = length(k);

% termo independente
b = zeros(nk,1);
b(1) = 1;

% matriz 'Akj'
A = zeros(nk,n);
for ik = 1:nk
    for ij = 0:n-1
        A(ik,ij+1) = ( ( 2*ij + 1 )/2 )^k(ik);
    end
end
A = 2*A;

% resolve o sistema
d = linsolve(A,b);

%-------------------------------------------%


% Display do vetor de coeficientes 'd'
fprintf('   1a derivada\n');
fprintf('(grid intercalado)\n');
fprintf('d=\n');% [');
for ik = 1:nk
    fprintf(' %14E\n',d(ik));
end
fprintf('--------------------\n');

%return


%% coeficientes para segunda derivada

%                                           %
% Monta e resolve o sistema 'Ax = b' -------%
%                                           %

% vetor k (0,2,...,2N)
k = 0:2:2*n;
nk = length(k);

% termo independente
b = zeros(nk,1);
b(2) = 2;

% matriz 'Akj'
A = zeros(nk,n);
%A(:,1) = ones(n,1);
A(1,1) = 1;
for ik = 1:nk
    for ij = 1:n%n-1
%         A(ik,ij+1) = 2*( ij + 1 )^k(ik);
        A(ik,ij+1) = 2*ij^k(ik);
    end
end
%A = 2*A;

% resolve o sistema
d = linsolve(A,b);

%-------------------------------------------%


% Display do vetor de coeficientes 'd'
fprintf('     2a derivada\n');
fprintf('(grid não intercalado)\n');
fprintf('d=\n');% [');
for ik = 1:nk
    fprintf(' %14E\n',d(ik));
end
%fprintf(' ]\n');


