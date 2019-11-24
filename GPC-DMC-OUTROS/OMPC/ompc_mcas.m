%%% Simulation of dual mode optimal predictive control
%%%
%%  [J,x,y,u,c,KSOMPC,F,t] = chap5_ompc_simulate_constraintsc(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist,umin,umax,Kxmax,xmax,rdmin,rdmax)
%%
%%   nc is the control horizon
%%   A, B, are the state space model parameters
%%   u = -Kx +c is nominal prediction for control law
%%
%%  Adds in constraint handling with constraints
%%  umin<u < umax   and    Kxmax*x < xmax
%%  rdmin< r-d < rdmax   Required to ensure feasibility
%%
%%   MCAS is F[x;c;r-d] <=t   or   N c + M x + V (r-d)<= t
%%
%%   xss = Kxr(r-d)    uss = Kur(r-d)
%%
%% Care must be taken to ensure the problem given is feasible

function [F,t,N,M,V,Kxr,Kur] = ompc_mcas(A,B,K,nc,umin,umax,Kxmax,xmax,rdmax,rdmin)


%%%% Estimate steady-state values
%%%  [xss;uss]=[Kxr;Kur](r-d)
M=inv([C,zeros(ny,nu);A-eye(nx),B]);
Kxr=M(1:nx,1:ny);
Kur=M(nx+1:nx+ny,1:ny);

%%%%% Define constraint matrices using invariant set methods on
%%%%%   Z= Psi Z  
%%%%%   u=-Kz Z  umin<u<umax   
%%%%%   Kxmax * x <xmax
%%%%%
%%%%% First define constraints at each sample as G*x<f
%%%%%
%%%%%  Find MAS as M x + N cfut + V*(r-d) <= f
Kz=[Kz,-K*Kxr-Kur];
G1=[-Kz;Kz;[Kxmax,zeros(size(Kxmax,1),nc*nu+nu)]];
f1=[umax;-umin;xmax]; 
G2=[zeros(2*nu,nx+nu*nc),[eye(nu);-eye(nu)]];
f2=[rdmax;-rdmin];
Phi=A-B*Kz(:,1:nx);
Psi=[Psi,[(eye(nx)-Phi)*Kxr;zeros(nc*nu,nu)]];
Psi=[Psi;[zeros(nu,nx+nc*nu),eye(nu)]];


[F,t]=findmas_tracking(Psi,G1,G2,f1,f2);
N=F(:,nx+1:nx+nu*nc);
M=F(:,1:nx);
V=F(:,nx+nu*nc+1:end);

