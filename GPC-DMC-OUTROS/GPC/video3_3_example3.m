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
ny=10;
nu=1;
lambda=1;
ny2=40;

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(a,b,ny);   
[H2,P2,Q2] = mpc_predmat(a,b,ny2);   
HH=tril(ones(nu,nu));

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);

%%% Control trajectory for a step change in r
Duopt = Pr;
yfut=H(:,1:nu)*Duopt;    %%% prediction within J
yfut2=H2(:,1:nu)*Duopt;  %%% Longer range prediction
uopt=HH*Duopt;

%%%% plotting
figure(1); clf reset
v=0:ny;v2=0:ny2;
aa=plot(v2,[0;yfut2],'b',[0 1 1 ny2],[0 0 1 1],'r--',...
    [ny,ny],[0,1.5],'g:',[0:nu-1,ny2],[uopt;uopt(end)],'m');
set(aa,'linewidth',2);
c=legend('Output','set point','horizon','input');
set(c,'fontsize',18)
title(['n_y =',num2str(ny),', n_u =',num2str(nu),', \lambda = ',num2str(lambda)],'fontsize',18)

%%%% plotting OL vs CL
%%% closed-loop
%%% Set point, disturbance and noise
ref = [zeros(1,5),ones(1,25)];
dist=[zeros(1,30)];
noise = [zeros(1,30)];
[y,u,Du,r] = mpc_simulate_noconstraints(b,a,nu,ny,lambda,1,ref,dist,noise);

figure(2); clf reset
v=0:ny;v2=0:ny2;
a=plot(v2,[0;yfut2],'b',[0 1 1 ny2],[0 0 1 1],'r--',...
    [ny,ny],[0,1.5],'g:',0:24,y(5:29),'b--');
xlim([0,25]);ylim([0 2])
set(a,'linewidth',2);
c=legend('open-loop','set point','horizon','closed-loop');
set(c,'fontsize',18)
title(['n_y =',num2str(ny),' n_u =',num2str(nu),', \lambda = ',num2str(lambda)],'fontsize',18)

