%%% Simulation of dual mode optimal predictive control using fixed
%%% constraint inequalities (no admissible sets)
%%%
%% 
%%%     x(k+1) = A x(k) + B u(k)            x0 is the initial condition
%%%     y(k) = C x(k) + D u(k) + dist       Note: Assumes D=0, dist unknown
%%
%% Uses observer to estimate the system states
%%
%%%%%   Constraints are
%%%%%   umin < u < umax    Kxmax * x <xmax
%%
%%%%%   The optimal cost is 
%%%%%     J = c'*SC*c + c'*SCX*(x-x_{ss}) + unconstrained optimal
%%
%%  Number of d.o.f.         nc
%% Horizon for constraint handling    npred
%%  Weighting matrices in J        Q, R
%%  Weighting matrices for terminal mode   Q2, R2
%%  reference trajectory     ref, r
%%  measurement noise        noise
%%  initial condition        x0
%%  
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

function [x,y,u,c,r] = sompc_simulate(A,B,C,D,Q,R,Q2,R2,nc,umin,umax,Kxmax,xmax,npred,x0,ref,dist,noise);

%%%% define all the matrices needed to implement constrained SOMPC
simdata=sompc_data(A,B,C,D,Q,R,Q2,R2,nc,umin,umax,Kxmax,xmax,npred);

L=simdata.L;K=simdata.K;Knew =simdata.Knew;
Ao=simdata.Ao;Bo=simdata.Bo;Co=simdata.Co;Do=simdata.Do;
CC=simdata.CC;dfixed=simdata.dfixed;dx0=simdata.dx0;
SX=simdata.SX;SXC=simdata.SXC;SC=simdata.SC;
Mx=simdata.Mx;Mu=simdata.Mu;Pr=simdata.Pr;

if norm(SXC)<1e-10; SXC=SXC*0;end
SC=(SC+SC')/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION

%%%%%%%%%% Initial Conditions 
nu=size(B,2);
c = zeros(nu*nc,2); u =zeros(nu,2);  x=[x0,x0];y=C*x;

z = [x;y*0];
r=ref;
runtime = size(ref,2)-1;

%%%%% Settings for quadratic program
opt = optimset('quadprog');
opt.Diagnostics='off';    %%%%% Switches of unwanted MATLAB displays
opt.LargeScale='off';     %%%%% However no warning of infeasibility
opt.Display='off';
opt.Algorithm='active-set';

Cost=0;


disp('********');
disp(['Number of inequalities is ',num2str(size(CC,1))]);
disp('***************');

for i=2:runtime;

%%%% Estiamte steady-state
distest=z(end-nu+1:end,i);
xss=Mx*(r(:,i)-distest);
uss=Mu*(r(:,i)-distest);

%%%%%%%%%%%%%%% CONSTRAINT HANDLING PART 
c(:,i+1) = c(:,i);  
[c(:,i),val,FLAG] =  quadprog(SC,[x(:,i)-xss]'*SXC/2,CC,dfixed+dx0*[z(:,i);r(:,i)],[],[],[],[],[],opt);
if FLAG==-2;disp('infeasible');
    
r=r*0;c=c*0;x=x*0;z=z*0;y=y*0; 
end
%%%%% Control law
u(:,i) = -Knew*z(:,i) + c(1:nu,i) + Pr*r(:,i);


%%%% Simulate model      
     x(:,i+1) = A*x(:,i) + B*u(:,i) ;
     y(:,i+1) = C*x(:,i+1) + dist(:,i+1);
%%%% Observer part  -removed for now
   %  z(:,i+1) = Ao*z(:,i) +Bo*u(:,i) + L*(y(:,i) + noise(:,i) - Co*z(:,i));
     z(:,i+1) = [x(:,i+1);dist(:,i+1)];
     
end



%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);
r = r(:,1:i+1);
noise = noise(:,1:i+1);
d = dist(:,1:i+1);

