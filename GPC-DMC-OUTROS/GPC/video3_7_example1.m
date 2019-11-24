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
ny=15:100;
nu=2;
lambda=1;


%%%% Generates a square H matrix for now


%%%% mpc_law computes entire proposed trajectory for several different nu
Nk=cell(1);Dk=Nk;Pr=Nk;

for ny=15:100;
    [H,P,Q] = mpc_predmat(a,b,ny);   
    [Nk{ny},Dk{ny},Pr{ny},S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);
    NN(ny)=Nk{ny}(1);
end

%%%% plotting
figure(1); clf reset

plot(NN,'b','linewidth',2);
title('Coefficients of Nk','fontsize',18)
xlabel('Prediction horizon');