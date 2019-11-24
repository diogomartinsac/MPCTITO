A =[0.9      0;
    0.2    0.1];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
Q=C'*C;
R=10;
nx=2;
nu=1;

%%% Horizon 1
nc=1;nuc=nu*nc;
[SX1,SC1,SXC1,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

%%% Horizon 2
nc=2;nuc=nu*nc;
[SX2,SC2,SXC2,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

%%% Horizon 3
nc=5;nuc=nu*nc;
[SX3,SC3,SXC3,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

[SX1,SX2,SX3]
SC1,SC2,SC3
[norm(SXC1),norm(SXC2),norm(SXC3)]


