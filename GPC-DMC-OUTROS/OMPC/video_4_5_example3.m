A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
D=zeros(2,2);
Q=C'*C;Q2=Q;
R=eye(2);R2=0.1*R;
nx=3;
nu=2;

%%% Horizon 1
nc=1;nuc=nu*nc;
[SX1,SC1,SXC1,Spsi,K1]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 2
nc=2;nuc=nu*nc;
[SX2,SC2,SXC2,Spsi,K2]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 3
nc=5;nuc=nu*nc;
[SX3,SC3,SXC3,Spsi,K3]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Horizon 4
nc=10;nuc=nu*nc;
[SX10,SC10,SXC10,Spsi,K10]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);

%%% Changes weights
nc=2;nuc=nu*nc;
[SX4,SC4,SXC4,Spsi]=chap4_suboptcost(A,B,C,D,Q,R,Q2,R2,nc);


[SX1,SX2,SX3,SX4]
SC1,SC2,SC3,SC4
[norm(SXC1),norm(SXC2),norm(SXC3),norm(SXC4)]

%%%% Control law u=-Ksompc x
KK=inv(SC1)*SXC1';
Ksompc=[K1+KK(1:nu,:)]
KK=inv(SC2)*SXC2';
Ksompc=[K2+KK(1:nu,:)]
KK=inv(SC3)*SXC3';
Ksompc=[K3+KK(1:nu,:)]
KK=inv(SC10)*SXC10';
Ksompc=[K10+KK(1:nu,:)]