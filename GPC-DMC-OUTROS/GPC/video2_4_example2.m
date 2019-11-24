%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -1.4 0.52 0.24];
b=[1,0.3,-0.2,1.1];
sizey=1;  %%% siso
ny=15;
nu=6;
lambda=3;

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(a,b,ny);   

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);

%%% use only first terms
Nk=Nk(1:sizey,:); Dk=Dk(1:sizey,:);
Pc = mpc_poles(a,b,Nk,Dk)

