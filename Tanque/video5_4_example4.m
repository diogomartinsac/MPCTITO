%SISO example of GPC simulation with constraints
%%  ay = bu is the model
a=[0 0.5927];% 0    0.5927
b=[1 -0.9759]; %1.0000   -0.9759
ny=10;  %% prediction horizon
nu=3;   %% input horizon
Wu=0.1; %%% weights on input
Wy=10;  %% weights on output
Dumax=0.1;  %% input rate limit
umax=9.5;   %% max input
umin=0.1;
ymax=1.2;
ymin=0;
sizeu=1; 
ref=[zeros(1,10),ones(1,50)];  %% target
dist=[zeros(1,20),0.1*ones(1,50)];   %%% output disturbance signal
noise=ref*0;  %% measurement noise

[y,u,Du,r] = mpc_simulate_outputconstraints(b,a,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);
