%MIMO example
A=[];B=[];
A(1,1:3:12) = poly([.5,.8,-.2]);
A(2,2:3:12) = poly([-.2,.9,.5]);
A(3,3:3:12) = poly([-.1,.4,.6]);
A(1,5:3:12) = [ -.2 .1 .02];
A(2,4:3:12) = [ .4 0 -.1];
A(3,4:3:12) = [0 .2 .2];
A(2,6:3:12) = [-1 .1 .3];
A(:,4:12) = A(:,4:12)*.8;
B = [.5 0.2 -.5 1 2 1;2 0 .3 -.8 .6 .5;0 .9 -.4 1 .3 .5];

ny=4;
nu=2;
sizeu=3;
[H,P,Q] = mpc_predmat(A,B,ny);
H=H(:,1:nu*sizeu);

Dumax=[1;0.4;0.6];
umax=[2.1;1.2;1.5];
umin=[-0.5;-0.4;-0.6];
ymax=[4;2;3];
ymin=[-1;-2;-2];
[CC,dd,du,ddu,dy]  = mpc_constraints2(Dumax,umax,umin,ymax,ymin,sizeu,nu,H,P,Q);
whos CC dd du ddu dy

