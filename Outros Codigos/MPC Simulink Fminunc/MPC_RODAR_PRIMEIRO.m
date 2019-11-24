clear all
clc
syms h1 h2 u1 u2 Q A1 A2 H c13 c23 c7 qout1 qout2 d1 d2
%Definição das entradas
u=[u1;u2];
%Definição dos estados
x=[h1;h2];        
%Equações diferenciais
%------------------------------------------
% Se h2 =< 0.405
qout1 = u1*sqrt(h1)/sqrt(c13*u1^2 + c23);
qout2 = c7*u2*sqrt(h2);

dh1dt = (Q - qout1)/A1;
dh2dt = (qout1 - qout2)/A2;  
sys1 = [dh1dt dh2dt];
%------------------------------------------
% Se h2 > 0.405
qout1 = u1*sqrt(h1-h2+H)/sqrt(c13*u1^2 + c23);
qout2 = c7*u2*sqrt(h2);

dh1dt = (Q - qout1)/A1;
dh2dt = (qout1 - qout2)/A2;  
sys2 = [dh1dt dh2dt];
%------------------------------------------
%% Linearização do Sistema - Matrizes com Simbólico
Am1=jacobian(sys1,x);
Bm1=jacobian(sys1,u);
Cm1= eye(2);
Dm1=zeros(size(Bm1));
%% -------------------------------------------
Am2=jacobian(sys2,x);
Bm2=jacobian(sys2,u);
Cm2= eye(2);
Dm2=zeros(size(Bm2));

%% Parâmetros
u1 = 0.2142;
u2 = 0.2411;
h1 = 0.405;
h2 = 0.405;
d1 = 0.13;
d2 = 0.06;
H = 0.405;
c13 = 3.4275e7;   
c23 = 0.9128e7;  
c7  = 2.7154e-4; 
A1 = (pi/4)*d1^4;
A2 = (pi/4)*d2^4;
Q = 150;

%Matrizes com valores numéricos

Am1=double(subs(Am1));
Bm1=double(subs(Bm1));
%-------------------------------------------
Am2=double(subs(Am2));
Bm2=double(subs(Bm2));

% State Space
Ts=0.1;
Dinamica_1 = ss(Am1,Bm1,Cm1,Dm1);
Dinamica_1d = c2d(Dinamica_1,Ts,'zoh');
%-------------------------------------------
Dinamica_2 = ss(Am2,Bm2,Cm2,Dm2);
Dinamica_2d = c2d(Dinamica_2,Ts,'zoh');

%MPC
nx = length(Dinamica_1d.A);    % Numero de estados
ny = size(Dinamica_1d.C,1);    % Numero de variaveis controladas
nu = size(Dinamica_1d.B,2);    % Numero de variaveis manipuladas
%
sysi1.A = [Dinamica_1d.A, Dinamica_1d.B; zeros(nu,nx), eye(nu)];
sysi1.B = [Dinamica_1d.B; eye(nu)]; sysi1.C = [Dinamica_1d.C, zeros(ny,nu)];
%--------------------------
sysi2.A = [Dinamica_2d.A, Dinamica_2d.B; zeros(nu,nx), eye(nu)];
sysi2.B = [Dinamica_2d.B; eye(nu)]; sysi2.C = [Dinamica_2d.C, zeros(ny,nu)];
nx = length(sysi1.A); % atualização do numero de estados
xk = zeros(nx,1);

% Parametros do controlador
Hp = 10; % horizonte de predicao
Hc = 3;  % Horizonte de controle

q = [1,1]; r = [0.0001/(u1^2),0.0001/(u2^2)];

% Matriz Phi
Phi1=[]; 
aux1=[]; 
for k=1:Hp
    Phi1=[Phi1; sysi1.C*sysi1.A^(k)];
    aux1=[aux1; sysi1.C*sysi1.A^(k-1)*sysi1.B];
end
% Matriz Theta
Theta1=aux1;
for k=1:Hc-1 
    aux1  = [zeros(ny,nu);aux1(1:(Hp-1)*ny,:)];
    Theta1 =[Theta1 aux1];
end
%-------------------------------------------
% Matriz Phi
Phi2=[]; 
aux2=[];
for k=1:Hp
    Phi2=[Phi2; sysi2.C*sysi2.A^(k)];
    aux2=[aux2; sysi2.C*sysi2.A^(k-1)*sysi2.B];
end
% Matriz Theta
Theta2=aux2;
for k=1:Hc-1 
    aux2  = [zeros(ny,nu);aux2(1:(Hp-1)*ny,:)];
    Theta2 =[Theta2 aux2];
end

% Calculo da Hessiana
Q1 = diag(repmat(q,1,Hp)); % Extensão da matriz de ponderacao
Q2 = diag(repmat(q,1,Hp)); % Extensão da matriz de ponderacao
R1 = diag(repmat(r,1,Hc)); % Extensão da matriz de ponderacao
R2 = diag(repmat(r,1,Hc)); % Extensão da matriz de ponderacao
H1 = Theta1'*Q1*Theta1 + R1;
H1 = (H1+H1')/2; % Ajuste de simetria numérica
%-------------------------------------------
H2 = Theta2'*Q2*Theta2 + R2;
H2 = (H2+H2')/2; % Ajuste de simetria numérica
%% Filtro de Kalman
Vd = 0.00001*eye(nx); % Variancia do modelo
Vn = 0.00001*eye(ny); % Variancia da medicao
xmk = xk; % Estado do modelo (valor inicial)
P = Vd;    % Variancia do erro (valor inicial)
P0 = 0.00001*eye(nx); 