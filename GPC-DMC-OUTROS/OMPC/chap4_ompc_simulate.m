%%% Simulation of dual mode optimal predictive control
%%%
%%  [J,x,y,u,c,KSOMPC] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime)
%%
%%   Q, R denote the weights in the actual cost function
%%   Q2, R2 are the weights used to find the terminal mode LQR feedback
%%   nc is the control horizon
%%   A, B,C,D are the state space model parameters
%%   x0 is the initial condition for the simulation
%%   J is the predicted cost at each sample
%%   c is the optimised perturbation at each sample
%%   x,y,u are states, outputs and inputs
%%   KSOMPC unconstrained feedback law

function [J,x,y,u,c,Ksompc] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime)

%%%%%%%%%% Initial Conditions 
nu=size(B,2);
c = zeros(nu*nc,2); u =zeros(nu,2);  
x=[x0,x0];
y=C*x;
runtime;
J=0;

%%%%%   The optimal predicted cost at any point 
%%%%%     J = c'*SC*c + 2*c'*SCX*x + x'*Sx*x
%%%% Control law parameters
[SX,SC,SXC,Spsi,K]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);
if norm(SXC)<1e-10; SXC=SXC*0;end
KK=inv(SC)*SXC';
Ksompc=[K+KK(1:nu,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION

for i=2:runtime;

%%%%% Control law
c(:,i) = KK*x(:,i);  
u(:,i) = -Ksompc*x(:,i);

%%%% Simulate model      
     x(:,i+1) = A*x(:,i) + B*u(:,i) ;
     y(:,i+1) = C*x(:,i+1);

%%% update cost
     J(i)=x(:,i)'*SX*x(:,i)+2*c(:,i)'*SXC'*x(:,i)+c(:,i)'*SC*c(:,i);
end

%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);
J(:,i+1)=J(:,i);
J(1)=J(2);


