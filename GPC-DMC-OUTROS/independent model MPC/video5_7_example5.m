%%% SISO Model and GPC parameters and T-filter
A=[1.4000   -0.1050   -0.1080; 2 0 0; 0 1 0];
B =[2; 0; 0]/10;
C =[ 0.5000    0.7500    0.5000]*10;
D=0;

%%% Tuning parameters
R=10;
ny=15;       %% Prediction horizon
nu=3;        %% Control horizon

Dumax=0.02;
umax=0.06;
umin=-0.05;

x0=[1;1;0]*0;  %% Initial condition
ref = [zeros(1,10),ones(1,50)];
dist= [zeros(1,20),ones(1,40)*0.6];
noise = [zeros(1,10),randn(1,50)*0];

%% figure 4-6
[x,y,u,r] = imgpc_simulate(A,B,C,D,R,ny,nu,umax,umin,Dumax,x0,ref,dist,noise);
%%
