%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -1.7 0.6];
b=[1 0.4];
sizey=1;  %%% siso
ny=20;
nu=5;
lambda=1;
ny2=60;

%%%% plotting OL vs CL
%%% closed-loop
%%% Set point, disturbance and noise
ref = [zeros(1,25),ones(1,ny2)];
dist=[zeros(1,ny2+25)];
noise = [zeros(1,ny2+25)];
[y,u,Du,r] = mpc_simulate_noconstraints_withff(b,a,nu,ny,lambda,1,ref,dist,noise);
[y2,u2,Du2,r] = mpc_simulate_noconstraints(b,a,nu,ny,lambda,1,ref,dist,noise);

figure(1); clf reset
v=0:ny;v2=0:ny2;
aa=plot([0 20 20 ny2],[0 0 1 1],'r--',...
    0:ny2-5,u(5:ny2),'g-',0:ny2-5,u2(5:ny2),'c-',0:ny2-5,y(5:ny2),'b-',0:ny2-5,y2(5:ny2),'m-');
xlim([0,40])
set(aa,'linewidth',2);
c=legend('set point','input with FF','input (no FF)','output with FF','output (no FF)');
set(c,'fontsize',14)
title(['n_y =',num2str(ny),' n_u =',num2str(nu),', \lambda = ',num2str(lambda)],'fontsize',18)

