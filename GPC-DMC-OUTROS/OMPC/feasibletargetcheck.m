%%%% Check that a desired target is feasible 
%%%%  for given input and state constraints
%%%%
%%%%%  gaps should all be negative

function gaps = feasibletargetcheck(A,B,C,umin,umax, Kxmax, xmax,rmax,rmin)

nx=size(A,1);
nu=size(B,2);
nK=size(Kxmax,1);

%%% determine expected steady-states
%% xss = Kxr*r  uss = Kur*r
Mat=inv([C,zeros(nu,nu);[A-eye(nx),B]]);
Kxr=Mat(1:nx,1:nu);
Kur=Mat(nx+1:end,1:nu);

%%%% Test for limits on target
%%%%  H*[rmax;rmin] <= h
%%%%  uses gap of epsilon =0.1 to ensure off  the boundary
H=[Kur,zeros(nu,nu);-Kur,zeros(nu,nu);zeros(nu,nu),Kur;zeros(nu,nu),-Kur;...
    Kxmax*Kxr,zeros(nK,nu);zeros(nK,nu),Kxmax*Kxr];
epsilon=0.01;
hh=[umax-epsilon;-umin+epsilon;umax-epsilon;-umin+epsilon;xmax-epsilon;xmax-epsilon];

gaps = H*[rmax;rmin]-hh;


