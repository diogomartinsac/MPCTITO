A =[0.9      0;
    0.2    0.1];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
Q=C'*C;Q2=Q;
R=eye(1); R2=15;
nx=2;
nu=1;

%%% Horizon 1
nc=1;nuc=nu*nc;
[SX1,SC1,SXC1,Spsi,K1]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 2
nc=2;nuc=nu*nc;
[SX2,SC2,SXC2,Spsi,K2]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 3
nc=3;nuc=nu*nc;
[SX3,SC3,SXC3,Spsi,K3]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 4
nc=10;nuc=nu*nc;
[SX10,SC10,SXC10,Spsi,K10]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

[SX1,SX2,SX3]
SC1,SC2,SC3
[norm(SXC1),norm(SXC2),norm(SXC3)]

%%%% Control law u=-Ksompc x
KK=inv(SC1)*SXC1';
Ksompc=[K1+KK(1,:)]
KK=inv(SC2)*SXC2';
Ksompc=[K2+KK(1,:)]
KK=inv(SC3)*SXC3';
Ksompc=[K3+KK(1,:)]
KK=inv(SC10)*SXC10';
Ksompc=[K10+KK(1:nu,:)]