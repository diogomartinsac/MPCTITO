function dhdt = conicalTank(t,h,wi_dot)
% Definitions of two-tank system at the Process Dynamics an Operantion
% Group at TU Dortmund.
% The system is two tanks connected by a pipe 5) and the available control
% inputs are the valves V3 and V7 that control the flow from tank B1 to
% tank B2 and the outflow from tank B2 respectively. A constant inflow
% through V10 of 150 l/h is assumed to for tank B1 and there is no
% additional inflow into tank B2.
% Parameters value of the two tank model
hM    =  2;                       % m
B       = 1.5;                     % m
C       = 4.0 ;                    % m
gamma = ((C/2) - (B/2))/hM;

rho    = 1000;         % kg/m^3
g       = 9.8;            % m/s^2
k       = 0.04 ;

c        = k*sqrt(rho*g);
% Planta

cplanta = 1.1*c;
% Computes the output flow of tank 1

beta = 1/(pi*gamma^2*(h+(B/2)/gamma)^2);

% ODES
dhdt = [beta*(wi_dot-cplanta*sqrt(h))];
end