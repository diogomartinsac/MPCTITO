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
ny=15;
nu=4;
ny2=25;
lambda=[0.01,0.1,1,10,100];

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(a,b,ny);   
[H2,P2,Q2] = mpc_predmat(a,b,ny2);   
HH=tril(ones(nu,nu));

%%%% mpc_law computes entire proposed trajectory for several different nu
Nk=cell(1);Dk=Nk;Pr=Nk;
Duopt=Nk;uopt=Nk;
yfut=Nk;yfut2=Nk;

for k=1:length(lambda);
    [Nk{k},Dk{k},Pr{k},S,X,Prlong] = mpc_law(H,P,Q,nu,lambda(k),1,sizey);
     HH=tril(ones(nu,nu));

%%% Control trajectory for a step change in r
Duopt{k} = Pr{k};
yfut{k}=H(:,1:nu)*Duopt{k};    %%% prediction within J
yfut2{k}=H2(:,1:nu)*Duopt{k};  %%% Longer range prediction
uopt{k}=HH*Duopt{k};
end

%%%% plotting
figure(1); clf reset
v=0:ny;v2=0:ny2;
plot([0,1,1,ny2],[0 0 1 1],'r--',[ny,ny],[0,1.5],'g:');hold on
for k=1:length(lambda);
    kk=(k-1)/(length(lambda)-1);
    plot(v2,[0;yfut2{k}],'linewidth',2,'color',[kk,0,1-kk]);
end
c=legend('set point','ny','\lambda=0.01','\lambda=0.1','\lambda=1','\lambda=10','\lambda=100');
set(c,'fontsize',18)
title(['Outputs with n_y =',num2str(ny),' nu = ',num2str(nu)],'fontsize',18)

%%%% plotting the inputs
figure(2); clf reset
v=0:ny;v2=0:ny2;
for k=1:length(lambda);
    kk=(k-1)/(length(lambda)-1);
    plot([0:nu-1,nu+2],[uopt{k};uopt{k}(end)],'linewidth',2,'color',[kk,0,1-kk]);
    hold on
end
c=legend('\lambda=0.01','\lambda=0.1','\lambda=1','\lambda=10','\lambda=100');
set(c,'fontsize',18)
title(['Inputs with n_y =',num2str(ny),' nu = ',num2str(nu)],'fontsize',18)
