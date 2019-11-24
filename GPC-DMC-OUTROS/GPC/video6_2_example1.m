%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -1.2 0.32];
b=[1,0.3];
sizey=1;  %%% siso
ny=20;
nu=20;
lambda=10;
ny2=30;
ref = [zeros(1,ny-1),ones(1,ny2-ny+1)];

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(a,b,ny);   
[H2,P2,Q2] = mpc_predmat(a,b,ny2);   
HH=tril(ones(nu,nu));

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);

%%% Control trajectory for a step change in r
Duopt = Prlong(:,end);
yfut=H(:,1:nu)*Duopt;    %%% prediction within J
yfut2=H2(:,1:nu)*Duopt;  %%% Longer range prediction
uopt=HH*Duopt;

%%%% plotting
figure(1); clf reset
v=0:ny;v2=0:ny2;
aa=plot(v2,[0;yfut2],'b',[0 ny ny ny2],[0 0 1 1],'r--',...
    [nu,nu],[-0.4,1.5],'g:',[0:nu-1,ny2],[uopt;uopt(end)],'m');
set(aa,'linewidth',2);
c=legend('Output','set point','n_u','input');
set(c,'fontsize',18)
title(['n_y =',num2str(ny),' n_u =',num2str(nu),' \lambda = ',num2str(lambda)],'fontsize',18)


