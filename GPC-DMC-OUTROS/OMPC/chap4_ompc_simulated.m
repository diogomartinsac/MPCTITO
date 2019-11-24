%%% Simulation of dual mode optimal predictive control with tracking
%%%
%%%  Uses an augmented model and observer to estimate states and 
%%%  disturbances and thus determine predictions
%%%
%%  [J,x,y,u,c,KSOMPC] = chap4_ompc_simulated(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise)
%%
%%   Q, R denote the weights in the actual cost function
%%   Q2, R2 are the weights used to find the terminal mode LQR feedback
%%   nc is the control horizon
%%   A, B,C, D are the independent state space model parameters
%%   Ap,Bp,Cp  represent the plant state space parameters
%%   x0 is the initial condition for the simulation
%%   J is the predicted cost at each sample
%%   c is the optimised perturbation at each sample
%%   x,y,u are states, outputs and inputs
%%   KSOMPC unconstrained feedback law
%%   ref a target signal
%%   output disturbance signal (assumes known for simplicity)

function [J,x,y,u,c,Ksompc,z] = chap4_ompc_simulatec(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise)

%%%%%%%%%% Initial Conditions 
nu=size(B,2);ny=size(C,1);nx=size(A,1);
c = zeros(nu*nc,2); u =zeros(nu,2);  
x=[x0,x0];z=x;xe=x;
y=C*x;
z=[x;y];
de=y;
runtime;
J=0;

%%%% Find observer and LQR feedback
%%%%   u = -K x - Kd d + Pr r   
%%%%   z=[x;d] is observer state
%%%%%  Equivalent to u-uss = -K(x-xss) is below
%%%%   u = -K x - Kd d + Pr r   
%%%%   z=[x;d] is observer state
[K,L,Ao,Bo,Co,Do,Kd,Pr,Mx,Mu] = ompc_observor(A,B,C,D,Q2,R2);

%%%%%   The optimal predicted cost at any point 
%%%%%     J = c'*SC*c + 2*c'*SCX*xhat + x'*Sx*xhat
%%%%%   xhat =x-xss
%%%% Control law parameters
[SX,SC,SXC,Spsi]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);
if norm(SXC)<1e-10; SXC=SXC*0;end
KK=inv(SC)*SXC';
Ksompc=[K+KK(1:nu,:)];


%%%% Estimate steady-state values and thus SOMPC feedback
%%%  [xss;uss]=[M1;M2](r-d)
%%%  u-uss = -Ksompc(x-xss)
%%%%   u = -Ksompc x - Kds d + Prs r
M=inv([C,zeros(ny,nu);A-eye(nx),B]);
M1=M(1:nx,1:ny);
M2=M(nx+1:nx+ny,1:ny);
Kds = Ksompc*M1+M2;
Prs=Kds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION

for i=2:runtime;

    
    %%%% collect state and disturbance estimates
    xe(:,i)=z(1:nx,i);
    de(:,i)=z(nx+1:end,i);
 
%%%%% Control law
u(:,i) = -Ksompc*xe(:,i) - Kds* de(:,i) + Prs*ref(:,i+1);  
     

%%%% Simulate process with disturbance      
     x(:,i+1) = Ap*x(:,i) + Bp*u(:,i) ;
     y(:,i+1) = Cp*x(:,i+1)+dist(:,i+1);
 
 %%%% Observer part
     z(:,i+1) = Ao*z(:,i) +Bo*u(:,i) + L*(y(:,i) - Co*z(:,i)+noise(:,i));

%%% update cost (based on deviation variables)
%%% To find current cost
xss(:,i)=M1*(ref(:,i+1)-de(:,i)); %% based on next target, current disturbance
uss(:,i)=M2*(ref(:,i+1)-de(:,i)); 
xhat(:,i)=x(:,i)-xss(:,i); 
c(:,i) = KK*xhat(:,i); 
J(i)=xhat(:,i)'*SX*xhat(:,i)+2*c(:,i)'*SXC'*xhat(:,i)+c(:,i)'*SC*c(:,i);

end

%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);
J(:,i+1)=J(:,i);
J(1)=J(2);


