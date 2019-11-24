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
% Acao integral
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


% q = [1.8*(1/h1^2) (1/h2^2)] ;
% r = [(1/u1^2) (1/u2^2)]/100;

% q = diag([50*(1/h1^2) 10*(1/h2^2) 60*(1/h1^2) 20*(1/h2^2)]);
% r = diag([(1/u1^2) (1/u2^2)]*5);

% q =[500*(1/0.405^2), 100*(1/0.405^2)] ;
% r =  [1*(1/0.2^2), 0.001*(1/0.2^2)]/100;
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
% Kf1=lqr(sysi1.A',sysi1.C',Vd,Vn)';
% sysKF1=ss((sysi1.A-Kf1*sysi1.C),[sysi1.B Kf1],eye(4),0*[sysi1.B Kf1]);
% %----------------
% Kf2=lqr(sysi2.A',sysi2.C',Vd,Vn)';
% sysKF2=ss((sysi2.A-Kf2*sysi2.C),[sysi2.B Kf2],eye(4),0*[sysi2.B Kf2]);
%% ARDUINO
% dlmwrite('invHESSIANA1.txt',H1^-1,'precision','%.8f');
% dlmwrite('fit_Q_theta1.txt',Phi1'*Q_*Theta1,'precision','%.8f');
% dlmwrite('Q_theta1.txt',Q_*Theta1,'precision','%.8f');
% %---------------------------------
% dlmwrite('invHESSIANA2.txt',H2^-1,'precision','%.8f');
% dlmwrite('fit_Q_theta2.txt',Phi2'*Q_*Theta2,'precision','%.8f');
% dlmwrite('Q_theta2.txt',Q_*Theta2,'precision','%.8f');
% 
% dlmwrite('AMODELO1.txt',sysi1.A,'precision','%.8f');
% dlmwrite('AMODELO2.txt',sysi2.A,'precision','%.8f');
% dlmwrite('BMODELO1.txt',sysi1.B,'precision','%.8f');
% dlmwrite('BMODELO2.txt',sysi2.B,'precision','%.8f');
% dlmwrite('CMODELO1.txt',sysi1.C,'precision','%.8f');
% dlmwrite('CMODELO2.txt',sysi2.C,'precision','%.8f');
% 
% type('invHESSIANA2.txt')
% type('Q_theta2.txt')
% type('fit_Q_theta2.txt')
% type('CMODELO2.txt')
% type('BMODELO2.txt')
% type('AMODELO2.txt')
% 
% type('invHESSIANA1.txt')
% type('Q_theta1.txt')
% type('fit_Q_theta1.txt')
% type('CMODELO1.txt')
% type('BMODELO1.txt')
% type('AMODELO1.txt')
% %% Restricoes
% umax = [0.95;0.95];  umin = [0.01; 0.01];
% Mtil=[];
% Itil=[];
% auxM=zeros(nu,Hc*nu);
% for in=1:Hc
%     auxM=[eye(nu) auxM(:,1:(Hc-1)*nu)];
%     Mtil=[Mtil;auxM];
%     Itil=[Itil;eye(nu)];
% end
% Ain = [Mtil;-Mtil];
% Bin = @(uk_1) [repmat(umax,Hc,1) - Itil*uk_1; Itil*uk_1 - repmat(umin,Hc,1)];
% 
% 
% du0 = zeros(nu,1);
%% Implementacao
nsim = 200;
% Modelo da planta
unc = 1; % incerteza de modelagem (unc = 1, indica que o modelo é perfeito)
planta.A = sysi1.A;
planta.B = sysi1.B*unc;
planta.C = sysi1.C*unc;
ysp=[0.2450;-0.1250];
xk = zeros(nx,1); uk_1 = zeros(nu,1);
for k = 1:nsim
 if k == 120
        ysp(1) = -0.1;
    end
    
    if k == 20
        ysp(2) = 0.2730;
    end    
 

    Ysp = repmat(ysp,Hp,1);
    
    % MPC
    %[du,J,ef] = fminunc(@(u) fob_in_mimo(xmk,u,Hp,Hc,sysi,ysp,q,r),du0);
    du = -H1^-1 * (xmk'*(Phi1'*Q1*Theta1) - Ysp'*Q1*Theta1)';
    % Insercao de disturbios externos
    if k==50
        de = 0;
    else
        de = 0;
    end
    
    % Planta                                  Variancia do modelo e medicao
    xk = planta.A*xmk + planta.B*(du(1:nu)+de);
    yk(:,k) = planta.C*xmk+ mvnrnd(zeros(ny,1),Vn)';
    uk(:,k) = uk_1 + du(1:nu);

    
    % Filtro de kalman - predicao
    xmk = (sysi1.A)*xmk + sysi1.B*(du(1:nu));
    ymk(:,k) = sysi1.C*xmk;
    P = sysi1.A*P*sysi1.A' + Vd;
    % Filtro de Kalman - correcao
    K = P*sysi1.C'/(sysi1.C*P*sysi1.C' + Vn);
    xmk = xmk + K*(yk(:,k) - ymk(:,k));
    P = (eye(nx) - K*sysi1.C)*P;
    
 % Atualizacao de variaveis
    du0 = [du(nu+1:end);zeros(nu,1)]; % Estimativa inicial do otmizador
    uk_1 = uk(:,k); % valor de u em k_1
    sp(:,k) = ysp;
end
%% Graficos
figure(1)
for iy = 1:ny
    subplot(ny,1,iy)
    plot(1:nsim,sp(iy,:),'--r',1:nsim,yk(iy,:),'.k',1:nsim,ymk(iy,:))
    ylabel(['y_' num2str(iy)])
    grid on
end
xlabel('k')
legend('Setpoint','Medicao','Modelo','Location','Best')
figure(2)
for iu = 1:nu
    subplot(nu,1,iu)
    stairs(1:nsim,uk(iu,:),'k')
    ylabel(['u_' num2str(iu)])
    grid on
end
xlabel('k')
