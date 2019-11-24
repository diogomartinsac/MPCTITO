%%% SISO Model and GPC parameters and T-filter
a=[1 -1.7 1.2];
b=[1 -1.4];
Tfilt=[1 -0.8];
sizey=1;  %%% siso
ny=12;
nu=4;
Wu=1;
Wy=1;
Dumax=0.2;
umax=1;
umin=-1.5;
ymax=1.1;
ymin=-0.1;
ref = [zeros(1,10),ones(1,50)];
dist=ref*0;
noise = [zeros(1,10),randn(1,50)*0.1];

%% figure 4-6
[y,u,Du,r] = mpc_simulate_tfiltoutputconstraints(b,a,Tfilt,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);
%% Figure 1-3
[y,u,Du,r] = mpc_simulate_outputconstraints(b,a,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);

