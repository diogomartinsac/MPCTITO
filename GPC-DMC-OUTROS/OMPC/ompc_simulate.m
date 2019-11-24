%%% Simulation of dual mode optimal predictive control
%%%
%%  [x,y,u,c,Sc,Sxc,Sx] = ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0)
%%
%%%%%   The optimal predicted cost at any point 
%%%%%     J = c'*SC*c + c'*SCX*x + x'*Sx*x
%%  

function [J,x,y,u,c] = ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime)

%%%%%%%%%% Initial Conditions 
nu=size(B,2);
c = zeros(nu*nc,2); u =zeros(nu,2);  
x=[x0,x0];
y=C*x;
runtime;
J=0;

%%%% Control law parameters
[SX,SC,SXC,Spsi,K]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

KK=inv(SC)*SXC';
Ksompc=[K1+KK(1:nu,:)];

if norm(SXC)<1e-10; SXC=SXC*0;end
SC=(SC+SC')/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION


for i=2:runtime;

%%%%% Control law
c(:,i) = KK*x(:,i);  
u(:,i) = -Ksompc*x(:,i);

%%%% Simulate model      
     x(:,i+1) = A*x(:,i) + B*u(:,i) ;
     y(:,i+1) = C*x(:,i+1) + dist(:,i+1);

%%% update cost
     J(i)=x(:,i)'*SX*x(:,i)+c(:,i)'*SXC*x(:,i)+c(:,i)'*SC*c(:,i);
end

%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);


