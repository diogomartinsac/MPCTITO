clc
clear all
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
sys1 = [dh1dt; dh2dt];
%------------------------------------------
% Se h2 > 0.405
qout1 = u1*sqrt(h1-h2+H)/sqrt(c13*u1^2 + c23);
qout2 = c7*u2*sqrt(h2);
dh1dt = (Q - qout1)/A1;
dh2dt = (qout1 - qout2)/A2;  
sys2 = [dh1dt; dh2dt];
%------------------------------------------
%% Linearização do Sistema - Matrizes com Simbólico
Am1=jacobian(sys1,x);
Bm1=jacobian(sys1,u);
Cm1= eye(2);
Dm1=0;
%% -------------------------------------------
Am2=jacobian(sys2,x);
Bm2=jacobian(sys2,u);
Cm2= eye(2);
Dm2=0;

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
Dinamica_1 = ss(Am1,Bm1,Cm1,Dm1);
%-------------------------------------------
Dinamica_2 = ss(Am2,Bm2,Cm2,Dm2);

%Ahat e Bhat ( Ação Integral )

Ahat1 = [Dinamica_1.A zeros(2,2)
        -Dinamica_1.C zeros(2,2)];
Bhat1 = [Dinamica_1.B
         zeros(size(Dinamica_1.C,1),size(Dinamica_1.B,2))];
%-------------------------------------------         
Ahat2 = [Dinamica_2.A zeros(2,2)
        -Dinamica_2.C zeros(2,2)];
Bhat2 = [Dinamica_2.B
         zeros(size(Dinamica_2.C,1),size(Dinamica_2.B,2))];

% LQR
% Q e R %
q = diag([50*(1/h1^2) 10*(1/h2^2) 60*(1/h1^2) 20*(1/h2^2)]);
r = diag([(1/u1^2) (1/u2^2)]*5);

% q = diag([500*1/h1^2, 100*1/h2^2 500*1/h1^2, 100*1/h2^2]) ;
% r = diag([1/u1^2, 1/u2^2]) ;
% q = diag([500*1/h1^2, 100*1/h2^2 600*1/h1^2, 200*1/h2^2]) ;
% r = diag([(1/u1^2), (1/u2^2)]) ;
% Ganhos
K1 = lqr(Ahat1,Bhat1,q,r);
Kp1 = K1(:,1:2);
Ki1 = K1(:,3:4);
%-------------------------------------------
K2 = lqr(Ahat2,Bhat2,q,r);
Kp2 = K2(:,1:2);
Ki2 = K2(:,3:4);
% %Kalman
% Vd=0.00001*eye(2)
% Vn=0.00001*eye(2)
% Kf=(lqr(Dinamica_1.A',Dinamica_1.C',Vd,Vn))'
% sysKF = ss((Dinamica_1.A-Kf*Dinamica_1.C),[Dinamica_1.B Kf],eye(2),0*[Dinamica_1.B Kf]);