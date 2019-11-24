%%% SISO Model and GPC parameters
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
ref = [zeros(1,10),ones(1,30)];
dist=ref*0;
noise = ref*0;
[y,u,Du,r] = mpc_simulate_tfilt(b,a,Tfilt,nu,ny,Wu,Wy,Dumax,umax,umin,ref,dist,noise);

