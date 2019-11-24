%%%%%  Typical data entry required for a 
%%%%%  SISO state space  example
%%%%%    x(k+1) = A x(k) + B u(k)
%%%%%    y(k) = C x(k) + D u(k) + dist     Note: Assumes D=0
%%%%%
%%%%%   Assumes J = sum (r-y)^2 +  (u(k+i-1)-uss) R (u(k+i-1)-uss)
%%%%%                uss the steady state input
%%%%%
%%%%%  Illustrates closed-loop simulations with constraint handling
%%%%%  
%%%%%
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% Model
A=[1.4000   -0.1050   -0.1080; 2 0 0; 0 1 0];
B =[2; 0; 0]/10;
C =[ 0.5000    0.7500    0.5000]*10;
D=0;

%%% Constraints
umax=.04;
umin=-.04;
Kxmax=[eye(3)];
xmax=[2;2;2];

%%% Tuning parameters
R=10;
Q=C'*C;
Q2=Q;R2=R;
npred=25;       %% Prediction horizon
nc=4;;        %% Control horizon

x0=[1;1;0]*0;  %% Initial condition

ref = [zeros(1,20),[1*ones(1,40)],ones(1,20)*2];    %%% Set point
dist = [zeros(1,35),-ones(1,45)*.5];   %%% disturbance
noise = [zeros(1,50),randn(1,30)*.04]; %%% noise

%%%%% Closed-loop simulation 
[x,y,u,c,r] = sompc_simulate(A,B,C,D,Q,R,Q2,R2,nc,umin,umax,Kxmax,xmax,npred,x0,ref,dist,noise);

figure(1); clf reset
v=2:length(y);runtime=length(y);
subplot(221);plot(v,y(:,v)','b','linewidth',2);hold on
plot(v,r(:,v)','k--');
title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(222);plot(v,u(:,v)','m','linewidth',2);hold on
plot([0,runtime],[umax,umax],'m--',[0,runtime],[umin,umin],'m--');
title('Inputs','fontsize',18)
subplot(223);plot(v,c(:,v)','linewidth',2);title('c_k for SOMPC','fontsize',18)
subplot(224);plot(v,dist(:,v)','linewidth',2);
title('Disturbance signal','fontsize',18)