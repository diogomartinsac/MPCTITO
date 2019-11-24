%%% Simulation of dual mode optimal predictive control
%%%
%%  [J,x,y,u,c,KSOMPC] = chap5_ompc_simulate_constraintsc(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist,umin,umax,Kxmax,xmax,rdmin,rdmax)
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
%%
%%   Assumes a time varying target signal and
%%   output disturbance signal 
%%
%%  Adds in constraint handling with constraints
%%  umin<u < umax   and    Kxmax*x < xmax
%%  rdmin< r-d < rdmax   Required to ensure feasibility
%%
%% Care must be taken to ensure the problem given is feasible

function [J,x,y,u,c,Ksompc] = chap5_ompc_simulate_constraintsc(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist,umin,umax,Kxmax,xmax,rdmax,rdmin)

%%%%%%%%%% Initial Conditions 
nu=size(B,2);
nx=size(A,1);
ny=size(C,1);
c = zeros(nu*nc,2); u =zeros(nu,2);  
x=[x0,x0];
y=C*x;
runtime;
J=0;

%%%%% The optimal predicted cost at any point 
%%%%%     J = c'*SC*c + 2*c'*SCX*x + x'*Sx*x
%%%%  Builds an autonomous model Z= Psi Z, u = -Kz Z  Z=[x;cfut];
%%%%
%%%% Control law parameters
[SX,SC,SXC,Spsi,K,Psi,Kz]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);
SC=(SC+SC')/2; %%% to avoid silly error messages
if norm(SXC)<1e-10; SXC=SXC*0;end
KK=inv(SC)*SXC';
Ksompc=[K+KK(1:nu,:)];

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
Psi=[Psi,[(eye(2)-Phi)*Kxr;zeros(nc*nu,nu)]];
Psi=[Psi;[zeros(nu,nx+nc*nu),eye(nu)]];

[F,t]=findmas_tracking(Psi,G1,G2,f1,f2);
N=F(:,nx+1:nx+nu*nc);
M=F(:,1:nx);
V=F(:,nx+nu*nc+1:end);

%%%%% Settings for quadratic program
opt = optimset('quadprog');
opt.Diagnostics='off';    %%%%% Switches of unwanted MATLAB displays
opt.LargeScale='off';     %%%%% However no warning of infeasibility
opt.Display='off';
opt.Algorithm='active-set';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION

for i=2:runtime;

rd(:,i)=ref(:,i+1)-dist(:,i);
xss(:,i)=Kxr*(ref(:,i+1)-dist(:,i)); %%  
uss(:,i)=Kur*(ref(:,i+1)-dist(:,i)); 
xhat(:,i)=x(:,i)-xss(:,i);

%%%%% Unconstrained control law
cfast(:,i) = KK*xhat(:,i);  
uhatfast(:,i) = -Ksompc*xhat(:,i);
ufast(:,i)=uhatfast(:,i)+uss(:,i);


%%%% constrained control law
%%%%  N c + M x +V(r-d) <=t
%%%%  J = c'*SC*c + 2*c'*SCX*xhat 
[cfut,vv,exitflag] = quadprog(SC,SXC'*xhat(:,i),N,t-M*x(:,i)-V*rd(:,i),[],[],[],[],[],opt);
if exitflag==-2;disp('No feasible solution');
    cfut=cfast(:,i);
%else; disp('looping')
    
end

c(:,i)=cfut;
uhat(:,i)=-K*xhat(:,i)+c(1:nu,i);
u(:,i)=uhat(:,i)+uss(:,i);


%%%% Simulate model      
     x(:,i+1) = A*x(:,i) + B*u(:,i) ;
     y(:,i+1) = C*x(:,i+1);

%%% update cost (based on deviation variables)
     J(i)=xhat(:,i)'*SX*xhat(:,i)+2*c(:,i)'*SXC'*xhat(:,i)+c(:,i)'*SC*c(:,i);

end

%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);
J(:,i+1)=J(:,i);
J(1)=J(2);


