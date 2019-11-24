%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% MIMO Model and GPC parameters
A=[];B=[];
A(1,1:3:12) = poly([.5,.8,-.2]);
A(2,2:3:12) = poly([-.2,.9,.5]);
A(3,3:3:12) = poly([-.1,.4,.6]);
A(1,5:3:12) = [ -.2 .1 .02];
A(2,4:3:12) = [ .4 0 -.1];
A(3,4:3:12) = [0 .2 .2];
A(2,6:3:12) = [-1 .1 .3];
A(:,4:12) = A(:,4:12)*.8;

B = [.5 0.2 -.5 1 2 1;2 0 .3 -.8 .6 .5;0 .9 -.4 1 .3 .5];
sizey=3;
ny=9;
nu=4;
lambda=0.2;

%%%% Generates a square H matrix for now
[H,P,Q] = mpc_predmat(A,B,ny);   

%%%% mpc_law computes entire proposed trajectory
[Nk,Dk,Pr,S,X,Prlong] = mpc_law(H,P,Q,nu,lambda,1,sizey);

%%% use only first terms
Nk=Nk(1:sizey,:); Dk=Dk(1:sizey,:);

