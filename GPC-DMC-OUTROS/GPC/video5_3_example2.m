%SISO example
a=[1 -0.8];
b=[0.4,0];
ny=5;
nu=3;
[H,P,Q] = mpc_predmat(a,b,ny);
H=H(:,1:nu);

Dumax=1;
umax=3;
umin=-2;
ymax=4;
ymin=-0.4;
sizeu=1; 
[CC,dd,du,ddu,dy]  = mpc_constraints2(Dumax,umax,umin,ymax,ymin,sizeu,nu,H,P,Q);
whos CC dd du ddu dy

