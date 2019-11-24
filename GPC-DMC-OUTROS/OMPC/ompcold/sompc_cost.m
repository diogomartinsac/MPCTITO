%%% Determine total predicted cost, including bits not depending on d.o.f.
%%% for a dual mode control law
%%% 
%%%% J = x SX x + x SXC c + c SC c
%%%%% Terminal mode control law and thus predictions based on Q, R
%%%%% Performance measured with Q2, R2
%%%     x(k+1) = A x(k) + B u(k)            
%%%     y(k) = C x(k) + D u(k) + dist       Note: Assumes D=0, dist unknown
%%%%%
%%%%% Regulation case only
%%%%%   
%%  
%%   [SX,SC,SXC,SSS]=sompc_cost(A,B,C,D,Q,R,Q2,R2,nc)
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)


function [SX,SC,SXC,SSS]=sompc_cost(A,B,C,D,Q,R,Q2,R2,nc)

nx = size(A,1);
nxc=nx*nc;
nu=size(B,2);
nuc=nu*nc;
%%% Split cost into 2 parts
%%  1. sum 1 to nc
%%  2. sum nc+1 to inf

%%%%%  Feedback loop is of the form  u = -Knew[x;d] + Pr*r  + c
%%%%%     Observor is   z = Ao*z +Bo*u + L*(y - Co*z);
[K,L,Ao,Bo,Co,Do,Knew,Pr] = ssmpc_observor(A,B,C,D,Q,R);
Phi=A-B*K;
%%%%%  x =  Pc1*c + Pz1*z + Pr1*r
%%%%%  u =  Pc2*c + Pz2*z + Pr2*r
%%%%%  y =  Pc3*c + Pz3*z + Pr3*r
[Pc1,Pc2,Pc3,Pz1,Pz2,Pz3,Pr1,Pr2,Pr3] = sompc_predclp(A,B,C,D,Knew,Pr,nc);
Px1=Pz1(1:nxc,1:nx);
Px2=Pz2(1:nuc,1:nx);
Pc1=Pc1(1:nxc,:);
Pc2=Pc2(1:nuc,:);

%%%%% Sum over transients 1 to nc
%%%%%  J1 = x Sx x + x Sxc c + c Sc c
for k=1:nc; v=(k-1)*nx+1:k*nx; v2=(k-1)*nu+1:k*nu;
    Qd(v,v)=Q2;
    Rd(v2,v2)=R2;
end

Sx = Px1'*Qd*Px1+Px2'*Rd*Px2;
Sxc = 2*(Px1'*Qd*Pc1+Px2'*Rd*Pc2);
Sc = Pc1'*Qd*Pc1+Pc2'*Rd*Pc2;

%% Sum over terminal mode [First capture x(nc) only]
Pxt=Px1(end-nx+1:end,:);
Pct=Pc1(end-nx+1:end,:);
Sxp = dlyap(Phi',Q2)-Q2;      %%% S = [Phi'*Q*Phi + Phi^2'*Q*Phi^2+...]
Sup = dlyap(Phi',K'*R2*K);   %%% S = [K'*R*K +Phi'K'*R*KPhi + Phi^2'K'*R*K Phi^2+...]
%%%% J2 = x(nc)[Sxp + Sup] x(nc) = x Sx2 x + x Sxc2 c + c Sc2 c
Sx2 = Pxt'*[Sxp+Sup]*Pxt;
Sxc2 = 2*Pxt'*[Sxp+Sup]*Pct;
Sc2 = Pct'*[Sxp+Sup]*Pct;

%%% Overall cost
%%%% J = x SX x + x SXC c + c SC c
SX=Sx+Sx2;
SXC=Sxc+Sxc2;
SC=Sc+Sc2;



