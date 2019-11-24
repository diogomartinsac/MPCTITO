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

%%% Horizon 1
nc=2;runtime=49;x0=[0;0];
[J,x,y,u,c,KSOMPC,z,ym] = chap4_ompc_simulatec(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);

figure(1); clf reset
v=2:length(y);
plot(v,y(v),'b','linewidth',2);title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,ref(v),'r','linewidth',2);
plot(v,dist(v),'g','linewidth',2);
plot(v,u(v),'m','linewidth',2);
plot(v,ym(:,v),'k--','linewidth',2)
legend('Output','target','disturbance','input','model');


