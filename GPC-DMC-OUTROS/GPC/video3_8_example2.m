%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -1.7 0.72];
b=[1 -2];
sizey=1;  %%% siso
ny=40;
nu=20;
lambda=1;
ny2=60;

%%%% plotting OL vs CL
%%% closed-loop
%%% Set point, disturbance and noise
ref = [zeros(1,5),ones(1,ny2)];
dist=[zeros(1,ny2+5)];
noise = [zeros(1,ny2+5)];
[y,u,Du,r] = mpc_simulate_noconstraints(b,a,nu,ny,lambda,1,ref,dist,noise);
[y2,u2,Du2,r2] = mpc_simulate_noconstraints(b,a,5,10,lambda,1,ref,dist,noise);

figure(2); clf reset
v=0:ny;v2=0:ny2;
aa=plot([0 1 1 ny2],[0 0 1 1],'r--',...
    0:ny2-5,u(5:ny2),'g-',0:ny2-5,y(5:ny2),'b-');
xlim([0,25])
set(aa,'linewidth',2);
c=legend('set point','input','output');
set(c,'fontsize',18)
title(['n_y =',num2str(ny),' n_u =',num2str(nu),', \lambda = ',num2str(lambda)],'fontsize',18)

figure(1); clf reset
v=0:ny;v2=0:ny2;
aa=plot([0 1 1 ny2],[0 0 1 1],'r--',...
    0:ny2-5,u(5:ny2),'g-',0:ny2-5,y(5:ny2),'b-',...
    0:ny2-5,u2(5:ny2),'g--o',0:ny2-5,y2(5:ny2),'b--o');

xlim([0,25])
set(aa,'linewidth',2);
c=legend('set point','input','output','input','output');
set(c,'fontsize',18)
title(['n_y =',num2str(10),' n_u =',num2str(5),', \lambda = ',num2str(lambda)],'fontsize',18)