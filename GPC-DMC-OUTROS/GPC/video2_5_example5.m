%%%%%  Typical data entry required for a 
%%%%%  SISO transfer function example
%%%%%    A(z) y  =  z^{-1}B(z) u
%%%%%
%%%%%  Illustrates overlaying closed-loop GPC simulations 
%%%%%  with NO constraint handling for different horizons
%%%%%  It should be obvious how to edit code to carry out any comparison
%%%%%  that the user may desire.
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
ny=[15,6,3,8,5];  % choices of prediction horizon
nu=[1 1 1 2 3];   % choices of input horizon

%%% Set point, disturbance and noise
ref = [zeros(1,5),ones(1,25)];
dist=[zeros(1,5),0*ones(1,25)];
noise = [zeros(1,15),randn(1,15)*0.0];

%%%%% Closed-loop simulation without and with a T-filter
Y=cell(1);U=Y;DU=Y;r=Y;
for k=1:length(ny);
[Y{k},U{k},DU{k},r{k}] = mpc_simulate_noconstraints(B,A,nu(k),ny(k),Wu,Wy,ref,dist,noise);
end

%%%%%% Overlaying the plots
figure(1); clf reset
nn=length(ny);
even=[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1]; %%% used for colours in plots

for k=1:nn;
time=0:size(U{k},2)-1;
Tun{k}=['Tuning (ny,nu,Wu)=(',num2str([ny(k),nu(k),Wu]),')'];

subplot(221);plot(time,Y{k},'Color',[1-k/nn,(k-1)/nn,even(k)]);hold on
xlabel(['GPC - Outputs and set-point']);
xlim([0,20]);ylim([0,1.2])
title('Overlaying GPC simulations for different tuning')
subplot(222);plot(time,DU{k}','Color',[1-k/nn,(k-1)/nn,even(k)]);hold on
xlabel(['GPC - Input increments ']);
xlim([0,20]);
subplot(223);plot(time,U{k}','Color',[1-k/nn,(k-1)/nn,even(k)]);hold on
xlabel(['GPC - Inputs']);
xlim([0,20]);ylim([0,0.5])

end
legend(Tun)
