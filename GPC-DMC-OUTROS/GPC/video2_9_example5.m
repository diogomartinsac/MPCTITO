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
A=[1 -1.2 0.32]; 
B=[1,.3];
Tfilt=[1 -1.6 0.64];
sizey=1;

%% Tuning parameters
Wu =1; % input weights
Wy=1;  % output weights
ny=15;  % prediction horizon
nu=3;   % input horizon

%%% Set point, disturbance and noise
ref = [zeros(1,5),0*ones(1,25)];
dist=[zeros(1,20),0*ones(1,10)];
noise = [zeros(1,5),randn(1,30)*0.02];

%%%%% Closed-loop simulation without and with a T-filter
[y,u,Du,r] = mpc_simulate_noconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise);
[yt,ut,Dut,r] = mpc_simulate_tfilt_noconstraints(B,A,Tfilt,nu,ny,Wu,Wy,ref,dist,noise);
    
%%%% Plotting
time=1:length(y);figure(1); clf reset
subplot(221);plot(time,y','b-',time,yt','r:',time,r','m--');ylim([-0.1,0.1])
xlabel(['Outputs and set-point']);
ss=legend('GPC','GPCT');set(ss,'fontsize',20);
subplot(222);plot(time,Du','b-',time,Dut','r:');ylim([-0.15,0.15])
xlabel(['Input increments']);
subplot(223);plot(time,u','b-',time,ut','r:');ylim([-0.05 0.05])
xlabel(['Inputs']);
subplot(224);plot(time,dist(time)','b',time,noise(time),'g');
xlabel(['Disturbance/noise']);

