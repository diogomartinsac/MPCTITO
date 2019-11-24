%%%% Predicting and control laws with a T-filter
%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -1.2 0.32];
b=[1,0.3];
Tfilt =[1 -0.8];
sizey=1;  %%% siso
ny=8;
nu=3;
lambda=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Find prediction matrices 
%%%%    yfut = H *Dufut + P*Dupast + Q*ypast
%%%%    yfut = H *Dufut + Pt*Dutpast + Qt*ytpast   %% filtered data
[H,P,Q] = mpc_predmat(a,b,ny);
[Pt,Qt] = mpc_predtfilt(H,P,Q,Tfilt,sizey,ny);

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);
[Nkt,Dkt,Prt,S,X,Prlongt] = mpc_law(H,Pt,Qt,nu,lambda,1,sizey);

%%% use only first terms
Nk=Nk(1:sizey,:); Dk=Dk(1:sizey,:);
Nkt=Nkt(1:sizey,:); Dkt=Dkt(1:sizey,:);

Pc = mpc_poles(a,b,Nk,Dk)
Pct = mpc_poles_tfilter(a,b,Nkt,Dkt,Tfilt)

roots(Pc),roots(Pct)