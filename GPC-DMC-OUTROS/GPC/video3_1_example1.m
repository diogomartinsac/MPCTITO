%%%%%  Typical data entry required for a 
%%%%%  SISO transfer function example
%%%%%    A(z) y  =  z^{-1}B(z) u
%%%%%
%%%%%  Illustrates closed-loop GPC simulations with NO constraint handling
%%%%%
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%% Model
A=[1 -1.2]; 
B=[1,.3]*.1;
sizey=1;

%% Tuning parameters
Wu =1; % input weights
Wy=1;  % output weights
ny=3;  % prediction horizon
nu=1;   % input horizon

%%% Set point, disturbance and noise
ref = [zeros(1,5),ones(1,15)];
dist=[zeros(1,30)];
noise = [zeros(1,30)];

%%%%% Closed-loop simulation without and with a T-filter
[y,u,Du,r] = mpc_simulate_noconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise);
    
%%%% plotting
figure(1); clf reset
v2=1:length(y);
aa=plot(v2,y,'b',v2,r,'r--',v2,u,'m');
set(aa,'linewidth',2);
c=legend('Output','set point','input');
set(c,'fontsize',18)
title(['n_y =',num2str(ny),', n_u =',num2str(nu),', \lambda = ',num2str(Wu)],'fontsize',18)


