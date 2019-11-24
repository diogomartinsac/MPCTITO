%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
%%% SISO Model and GPC parameters
a=[1 -1.7 0.6];
b=[1 0.4];
sizey=1;  %%% siso
ny=15;
nu2=5;
lambda=1;
ny2=25;

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(a,b,ny);   
[H2,P2,Q2] = mpc_predmat(a,b,ny2);   
HH=tril(ones(nu,nu));

%%%% mpc_law computes entire proposed trajectory for several different nu
Nk=cell(1);Dk=Nk;Pr=Nk;
Duopt=Nk;uopt=Nk;
yfut=Nk;yfut2=Nk;

for nu=1:nu2;
    [Nk{nu},Dk{nu},Pr{nu},S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);
     HH=tril(ones(nu,nu));

%%% Control trajectory for a step change in r
Duopt{nu} = Pr{nu};
yfut{nu}=H(:,1:nu)*Duopt{nu};    %%% prediction within J
yfut2{nu}=H2(:,1:nu)*Duopt{nu};  %%% Longer range prediction
uopt{nu}=HH*Duopt{nu};
end

%%%% plotting
figure(1); clf reset
v=0:ny;v2=0:ny2;
plot([0,1,1,ny2],[0 0 1 1],'r--',[ny,ny],[0,1.5],'g:');hold on
for nu=1:nu2;
    k=(nu-1)/(nu2-1);
    plot(v2,[0;yfut2{nu}],'linewidth',2,'color',[k,0,1-k]);
end
c=legend('set point','ny','nu=1','nu=2','nu=3','nu=4','nu=5');
set(c,'fontsize',18)
title(['Predictions with n_y =',num2str(ny),' \lambda = ',num2str(lambda)],'fontsize',18)

%%%% plotting the inputs
figure(2); clf reset
v=0:ny;v2=0:ny2;
plot([nu2-1,nu2-1],[-0.2,0.4],'g:');hold on
for nu=1:nu2;
    k=(nu-1)/(nu2-1);
    plot([0:nu-1,nu2+2],[uopt{nu};uopt{nu}(end)],'linewidth',2,'color',[k,0,1-k]);
end
c=legend('nu(max)','nu=1','nu=2','nu=3','nu=4','nu=5');
set(c,'fontsize',18)
title(['Inputs with n_y =',num2str(ny),' \lambda = ',num2str(lambda)],'fontsize',18)
