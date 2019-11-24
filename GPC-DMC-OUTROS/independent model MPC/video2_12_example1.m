%%%%%  Typical data entry required for a 
%%%%%  SISO transfer function example
%%%%%    A(z) y  =  z^{-1}B(z) u
%%%%%
%%%%%  Illustrates closed-loop GPC simulations with NO constraint handling
%%%%%  and using independent models for prediction
%%%%%
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%% Model
A=[1 -1.2 0.32]; 
B=[1,.3];
sizey=1;

%% Tuning parameters
Wu =1; % input weights
Wy=1;  % output weights
ny=15;  % prediction horizon
nu=3;   % input horizon

%%%%% Closed-loop simulation 1
ref = [zeros(1,5),ones(1,25)];
dist=[zeros(1,15),0.2*ones(1,15)];
noise = [zeros(1,30)];
[y,u,Du,r] = imgpc_tf_simulate_noconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison with GPC
%%%%% Closed-loop simulation 2
ref = [zeros(1,5),0*ones(1,25)];
dist=[zeros(1,30)];
noise = [zeros(1,10),randn(1,20)*0.01];
[yim,uim,Duim,r] = imgpc_tf_simulate_noconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise);
%%%% NOTE FILE BELOW IS IN A DIFFERENT FOLDER
[y,u,Du,r] = mpc_simulate_noconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise);

figure(1);clf reset
t=1:length(u);
subplot(211); plot(t,u,'b',t,uim,'r');
title('Inputs','fontsize',20)
legend('GPC','IMPGPC');
subplot(212);
plot(t,Du,'b',t,Duim,'r');
title('Input increments','fontsize',20);
