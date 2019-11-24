%%%%%   THIS IS A FUNCTION FILE TO DETERMINE HOW SENSITIVITY 
%%%%%   FOR GPC CONTROL LAWS
%%%%     COMPARES WITH AND WITHOUT A T-filter.
%%%%    ONLY COMPUTES SENSITIVITY TO MEASUREMENT NOISE (INPUT and OUTPUT)
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)
%%
%% [S,ST,w,Sy,SyT]=mpc_sensitivity(a,b,Tfilt,Wu,Wy,ny,nu,sizey)
%%
%% Model is ay=bu, with Tfilter given as Tfilt
%% Wu, Wy are eights in J 
%% ny,nu GPC are horizons
%%
%% ONLY WORKS FOR SISO CASES

function [mag1,mag2,w,mag1y,mag2y]=mpc_sensitivity(a,b,Tfilt,lambda,Wy,ny,nu)

sizey=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Find prediction matrices 
%%%%    yfut = H *Dufut + P*Dupast + Q*ypast
%%%%    yfut = H *Dufut + Pt*Dutpast + Qt*ytpast   %% filtered data
[H,P,Q] = mpc_predmat(a,b,ny);
[Pt,Qt] = mpc_predtfilt(H,P,Q,Tfilt,sizey,ny);

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,Wy,sizey);
[Nkt,Dkt,Prt,S,X,Prlongt] = mpc_law(H,Pt,Qt,nu,lambda,Wy,sizey);

%%% use only first terms
Nk=Nk(1:sizey,:); Dk=Dk(1:sizey,:);
Nkt=Nkt(1:sizey,:); Dkt=Dkt(1:sizey,:);

Pc = mpc_poles(a,b,Nk,Dk);
[Pct,TD] = mpc_poles_tfilter(a,b,Nkt,Dkt,Tfilt)

%%%% not being careful about exact value of delay as want norm only
S=tf(conv(a,Nk),Pc,1);
ST=tf(conv(a,Nkt),conv(Pc,Tfilt),1);
Sy=tf(conv(conv(a,Dk),[1 -1]),Pc,1);
SyT=tf(conv(conv(a,TD),[1 -1]),conv(Pc,Tfilt),1);

%%%% Calculate gains on same spectrum
[mag1,ph,w]=bode(S);
[mag2]=bode(ST,w);  %% with T-filter
[mag1y]=bode(Sy,w);
[mag2y]=bode(SyT,w); %% with T-filter


