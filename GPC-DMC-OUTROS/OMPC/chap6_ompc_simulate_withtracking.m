%%% Simulation of dual mode optimal predictive control with tracking AND
%%% ADVANCE KNOWLEDGE OF SET POINT
%%%
%%%  Uses an augmented model and observer to estimate states and 
%%%  disturbances and thus determine predictions
%%%
%%  [x,y,u,c,z] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise)
%%
%%   Q, R denote the weights in the actual cost function
%%   Q2, R2 are the weights used to find the terminal mode LQR feedback
%%   nc is the control horizon
%%   na is the advance knowledge of the set point
%%   A, B,C, D are the independent state space model parameters
%%   Ap,Bp,Cp  represent the plant state space parameters
%%   x0 is the initial condition for the simulation
%%   J is the predicted cost at each sample
%%   c is the optimised perturbation at each sample
%%   x,y,u are states, outputs and inputs
%%   ref a target signal
%%   output disturbance signal (assumes known for simplicity)

function [x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise)

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

%%%% Estimate steady-state values and thus terminal feedback without c
%%%  [xss;uss]=[M1;M2](r-d)
%%%  u-uss = -K(x-xss)
%%%%   u = -K x +Kds(r-d)
M=inv([C,zeros(ny,nu);A-eye(nx),B]);
Kxr=M(1:nx,1:ny);
Kur=M(nx+1:nx+ny,1:ny);


%%%%%   The optimal predicted cost at any point 
%%%% J = x SX x+2*xSXC c+c SC c+2*x*Sxr*rfut+2*c*Scr*rfut+rfut*Sr*rfut
%%%%%
%%%%%      with control law is given as u-uss = Kzss[x;c;r-d]
%%%
%%%%   Optimal perturbation c is given from
%%%%   cfut = -inv(SC)[SXC'*x+Scr*(rfut-d)]
%%%%   cfut = -Ksompc*x -Pr*(rfut-d)
[SX,SC,SXC,Scr,Kzss,Spsi]=chap6_ompccost_withtracking(A,B,C,D,Q,R,Q2,R2,nc,na);
if norm(SXC)<1e-10; SXC=SXC*0;end
Ksompc=inv(SC)*SXC';
Pr=inv(SC)*Scr;


%%% define matrix for creating [dist;dist;dist,...];
LL=[];
for k=1:na;
    LL=[LL;eye(nu)];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   SIMULATION

for i=2:runtime;

    
    %%%% collect state and disturbance estimates
    xe(:,i)=z(1:nx,i);
    de(:,i)=z(nx+1:end,i);
    uss(:,i)=Kur*(ref(:,i+1)-de(:,i)); 

    %%%% Vector of future targets (rfut - dist)
    rfut=ref(:,i+1:i+na);
    rfut_dist=rfut(:)-LL*de(:,i);
    
     %%% Optimal perturbations
     cfut=-Ksompc*xe(:,i)-Pr*rfut_dist;
     
%%%%% Control law
%%% u-uss = Kzss*[x;cfut;rfut_dist] 

u(:,i) = Kzss*[xe(:,i);cfut;rfut_dist]+uss(:,i);  
     

%%%% Simulate process with disturbance      
     x(:,i+1) = Ap*x(:,i) + Bp*u(:,i) ;
     y(:,i+1) = Cp*x(:,i+1)+dist(:,i+1);
 
 %%%% Observer part
     z(:,i+1) = Ao*z(:,i) +Bo*u(:,i) + L*(y(:,i) - Co*z(:,i)+noise(:,i));


xss(:,i)=Kxr*(ref(:,i+1)-de(:,i)); %% based on next target, current disturbance
c(:,i) = cfut; 

end

%%%% Ensure all variables have conformal lengths
u(:,i+1) = u(:,i);  
c(:,i+1)=c(:,i);


