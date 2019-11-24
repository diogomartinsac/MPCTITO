%%%  Determine suboptimal MPC cost, including bits not depending on d.o.f.
%%%  nc is the control horizon.
%%%  na is the advance knowledge of the target trajectory
%%%
%%%  Assumes terminal mode not matched to the LQR regulator.
%%%
%%%  Performance index based on    J =sum x'Q x+u'R u
%%% Terminal Control law u-uss=-K(x-xss) and predictions based on Q2, R2
%%%
%%%  Assume model x(k+1)=Ax(k)+Bu(k), y= Cx+Du
%%%
%%%  Overall cost
%%%% J = x SX x+2*xSXC c+c SC c+2*x*Sxr*rfut+2*c*Scr*rfut+rfut*Sr*rfut
%%%%  with control law is given as u-uss = Kzss[x;c;r-d]
%%%%
%%%%   [SX,SC,SXC,Scr,Kzss,Spsi,K,Psi,Kz]=chap6_ompccost_withtracking(A,B,C,D,Q,R,Q2,R2,nc,na)
%%%%
%%%%  Builds an autonomous model Z= Psi Z, u = Kz Z  Z=[x;cfut;rfut-d];
%%%%  Uses to form a cost based on deviations     x-xss and u-uss.
%%%%
%%%% Code by JA Rossiter (2014)

function [SX,SC,SXC,Scr,Kzss,Spsi,K,Psi,Kz]=chap6_ompccost_withtracking(A,B,C,D,Q,R,Q2,R2,nc,na)

nx = size(A,1);
nxc=nx*nc;
nu=size(B,2);
ny=nu;
nuc=nu*nc;

%%%%%  Feedback loop is of the form  u = -Kx+c
%%%%%  Find terminal mode feedback (uses LQR based on Q2, R2) 
[K] = dlqr(A,B,Q2,R2);
Phi=A-B*K;


%%%% Estimate steady-state values and thus terminal feedback without c
%%%  [xss;uss]=[Kxr;Kur](r-d)
%%%  u-uss = -K(x-xss)
M=inv([C,zeros(ny,nu);A-eye(nx),B]);
Kxr=M(1:nx,1:ny);
Kur=M(nx+1:nx+ny,1:ny);

%%% Build autonomous model without tracking
%%%%   Z= Psi Z, u = -Kz Z  Z=[x;cfut];

ID=diag(ones(1,(nc-1)*nu));
ID=[zeros((nc-1)*nu,nu),ID];
ID=[ID;zeros(nu,nuc)];
Psi=[A-B*K,[B,zeros(nx,nuc-nu)];zeros(nuc,nx),ID];
Kz = [K,-eye(nu),zeros(nu,nuc-nu)];

%%% Build autonomous model with tracking
%%%%  Z= Psi Z, u = Kz Z  Z=[x;cfut;rfut-d];

DR=diag(ones(1,(na-1)*nu));
DR=[zeros((na-1)*nu,nu),DR];
DR=[DR;zeros(nu,nu*na)];
DR((na-1)*nu+1:na*nu,(na-1)*nu+1:na*nu)=eye(nu);
Psi=[Psi,[[(eye(nx)-Phi)*Kxr,zeros(nx,(na-1)*nu)];zeros(nc*nu,na*nu)]];
Psi=[Psi;[zeros(na*nu,nx+nc*nu),DR]];

%%% Express input in terms of augmented state
%%%  u-uss=-K(x-xss) = -Kx +K xss +uss
%%%     u = [-K,0,K*Kxr+Kur][x;c;r-d]
%%%  u = -Kz*[x;c;r-d]
Kz=[Kz,-K*Kxr-Kur,zeros(nu,(na-1)*nu)];

%%%% Substitute the autonomous model into the performance index
%%%% Solve for the cost parameters using lyapunov
Gamma=[eye(nx),zeros(nx,nu*nc),-Kxr,zeros(nx,(na-1)*nu)];
Kzss=-Kz+[zeros(nu,nx+nu*nc),-Kur,zeros(nu,(na-1)*nu)];
W=Psi'*Gamma'*Q*Gamma*Psi+Kzss'*R*Kzss;
%%%% NOTE: dlyap needs Psi to be strictly stable, so a work around is to set
%%%% the terminal target dynamics to change very very slowly as this has minimal
%%%% impact on the performance index
Psi(end-nu+1:end,end-nu+1:end)=eye(nu)*0.99999999;
Spsi=dlyap(Psi',W);

%%% Split cost into parts
%%%% J = [x,c,r] Spsi [x;c,r]
%%%% J = x SX x+2*xSXC c+c SC c+2*x*Sxr*rfut+2*c*Scr*rfut+rfut*Sr*rfut
SX=Spsi(1:nx,1:nx);
SXC=Spsi(1:nx,nx+1:nx+nu*nc);
SC=Spsi(nx+1:nx+nc*nu,nx+1:nx+nc*nu);
Scr=Spsi(nx+1:nx+nu*nc,nx+nu*nc+1:end);





