A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
D=zeros(2,2);
Q=C'*C;
R=eye(2);
nx=3;
nu=2;

%%% Horizon 1
nc=1;nuc=nu*nc;
[SX1,SC1,SXC1,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

%%% Horizon 2
nc=2;nuc=nu*nc;
[SX2,SC2,SXC2,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

%%% Horizon 3
nc=5;nuc=nu*nc;
[SX3,SC3,SXC3,Spsi]=chap4_cost(A,B,C,D,Q,R,nc);

%%% Changes weights
nc=2;nuc=nu*nc;
[SX4,SC4,SXC4,Spsi]=chap4_cost(A,B,C,D,Q,R*10,nc);


[SX1,SX2,SX3,SX4]
SC1,SC2,SC3,SC4
[norm(SXC1),norm(SXC2),norm(SXC3),norm(SXC4)]


