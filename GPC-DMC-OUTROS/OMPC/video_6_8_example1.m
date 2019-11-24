A =[0.8      0;
    0.2    0.2];
B =[0.2;0.8];
C =[1.9 -1];
Ap=A;Bp=B;Cp=C;
Ap(2,1)=0.3;Ap(2,2)=0.3;Bp(1,1)=0.4;Cp(1,1)=2;
D=0;
Q=C'*C;Q2=Q;
R=eye(1); R2=R;
nx=2;
nu=1;
ref=[zeros(1,10),ones(1,50)];
dist=[zeros(1,30),-0.5*ones(1,30)];
noise=dist*0;
nc=5;
na=5;
runtime=9;x0=[0;0];
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL

na=3;
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL

na=1;
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL
