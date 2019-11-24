A =[0.8      0;
    0.2    0.2];
B =[0.2;0.8];
C =[1.9 -1];
Ap=A;Bp=B;Cp=C;
%Ap(2,1)=0.3;Ap(2,2)=0.3;Bp(1,1)=0.4;Cp(1,1)=2;
D=0;
Q=C'*C;Q2=Q;
R=eye(1); R2=0.1*R;
nx=2;
nu=1;
ref=[zeros(1,10),ones(1,80)];
dist=[zeros(1,30),-0.5*ones(1,60)];
noise=[zeros(1,60),randn(1,30)*0.1];

%%% Horizon 1
nc=2;runtime=79;x0=[0;0];
[J,x,y,u,c] = chap4_ompc_simulatec(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
[J2,x2,y2,u2,c2] = chap4_ompc_simulated(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);

figure(1); clf reset
v=2:length(y);
subplot(211)
plot(v,y(v),'b','linewidth',2);title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,y2(v),'r','linewidth',2);
title('Outputs')
subplot(212)
plot(v,u(v),'b','linewidth',2);hold on
plot(v,u2(:,v),'r','linewidth',2)
legend('Independent','Observer');
title('inputs')



