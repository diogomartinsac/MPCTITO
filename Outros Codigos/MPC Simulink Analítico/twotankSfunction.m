function [sys,x0] = twotankSfunction(t,x,u,flag)
% Definitions of two-tank system at the Process Dynamics an Operantion
% Group at TU Dortmund.
% The system is two tanks connected by a pipe 5) and the available control
% inputs are the valves V3 and V7 that control the flow from tank B1 to
% tank B2 and the outflow from tank B2 respectively. A constant inflow
% through V10 of 150 l/h is assumed to for tank B1 and there is no
% additional inflow into tank B2.
% Parameters value of the two tank model
d = [0.13; 0.06]; % m - Tank 1 and 2 diameter
H = 0.405;        % m - Conection heigth
c13 = 3.4275e7;   % s2/m5 - valve constant
c23 = 0.9128e7;   % s2/m5 - valve constant
c7  = 2.7154e-4;  % s2/m5 - valve constant
A = (pi/4)*d.^4;  % m^2 tank1 and tank2 area

if abs(flag) == 1
   x(x <= 0) = 0; x(x >= 1) = 1;
   % Cálculo das derivadas dos estados
   Q = u(3)/(3.6e6); h = x;
    % Computes the output flow of tank 1
    if h(2) <= H
        q_out(1) = u(1)*sqrt(h(1))/sqrt(c13*u(1)^2 + c23);
    else
        q_out(1) = u(1)*sqrt(h(1)-h(2)+H)/sqrt(c13*u(1)^2 + c23);
    end
    % Computes the output flow of tank 2
    q_out(2) = c7*u(2)*sqrt(h(2));

    % ODES
    dhdt = [(Q - q_out(1))/A(1)
            (q_out(1) - q_out(2))/A(2)];
        
   sys = dhdt       ;

elseif abs(flag) == 3
   x(x <= 0) = 0; x(x >= 1) = 1;
   sys = x;

elseif flag == 0
% Inicialização do sistema
  x0       = [0.4050, 0.4050];
  sys(1,1) =  2 ; 
  sys(2,1) =  0 ; 
  sys(3,1) =  2 ; 
  sys(4,1) =  3 ; 
  sys(5,1) =  0 ; 
  sys(6,1) =  1 ; 

else
  sys = [];  
end

% Autor:  Ricardo Kalid
%         kalid@ufba.br

% Fim deste arquivo
