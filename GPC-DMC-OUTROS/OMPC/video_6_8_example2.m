A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
Ap=A;Bp=B;Cp=C;
D=zeros(2,2);
Q=C'*C;Q2=Q;
R=eye(2);R2=R;
nx=3;
nu=2;
na=5;
nc=5;
ref=[zeros(2,8),[ones(1,100);zeros(1,10),0.5*ones(1,90)]];
dist=[zeros(2,70),ones(2,35)*0.0];
noise=dist*0;

runtime=9;x0=[0;0;0];
na=5;
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL

na=3;
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL

na=1;
[x,y,u,c,z,Pr,LL] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
Pr,LL,Pr*LL
