%SISO example of GPC simulation with constraints
%%  ay = bu is the model
a=[1 -0.8];
b=[0.1,0];
Tfilt=[1 -0.8];
ny=10;  %% prediction horizon
nu=2;   %% input horizon
Wu=1; %%% weights on input
Wy=1;  %% weights on output
Dumax=0.5;  %% input rate limit
umax=3;   %% max input
umin=-2;
ymax=1.2;
ymin=-0.2;
sizeu=1; 
ref=[zeros(1,10),ones(1,40)];  %% target
dist=[zeros(1,20),0.4*ones(1,30)];   %%% output disturbance signal
noise=ref*0;  %% measurement noise

[y,u,Du,r] = mpc_simulate_tfiltoutputconstraints(b,a,Tfilt,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);

[y,u,Du,r] = imgpc_tf_simulate_constraints(b,a,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);