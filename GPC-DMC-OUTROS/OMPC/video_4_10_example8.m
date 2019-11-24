A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
Ap=A;Bp=B;Cp=C;
Ap(1,1)=0.9;Ap(2,2)=0.12;Ap(2,3)=0.1;Bp(1,2)=-0.1;Bp(2,2)=0.15;Cp(1,2)=12;
D=zeros(2,2);
Q=C'*C;Q2=Q;
R=eye(2);R2=0.1*R;
  %R2=R; %to get OMPC
nx=3;
nu=2;
ref=[zeros(2,5),[ones(1,130);zeros(1,30),0.5*ones(1,100)]];
dist=[zeros(2,70),ones(2,65)*0.4];
noise=[zeros(2,100),randn(2,30)*0.1];

%%% Horizon 1
nc=5;runtime=129;x0=[0;0;0];
[J,x,y,u,c] = chap4_ompc_simulatec(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
[J2,x2,y2,u2,c2] = chap4_ompc_simulated(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);

%%%%% plotting
figure(1); clf reset
v=2:size(y,2);
subplot(211)
plot(v,y(:,v)','b','linewidth',2);title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,y2(:,v)','r','linewidth',2);
title('Outputs')
subplot(212)
plot(v,u(:,v)','b','linewidth',2);hold on
plot(v,u2(:,v)','r','linewidth',2)
legend('IM 1','IM 2','Observer');
title('Inputs')



