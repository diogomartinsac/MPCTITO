A =[0.8      0;
    0.2    0.2];
B =[0.2;0.8];
C =[1.9 -1];
D=0;
Q=C'*C;Q2=Q;
R=0.1*eye(1); R2=R;
nx=2;
nu=1;
ref=[zeros(1,10),ones(1,50)];
dist=[zeros(1,30),-0.5*ones(1,30)];

%%% Horizon 1
nc=2;runtime=49;x0=[0.2;-0.4];
[J,x,y,u,c,KSOMPC] = chap4_ompc_simulateb(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist);

figure(1); clf reset
v=2:length(y);
subplot(311);plot(v,y(v),'b','linewidth',2);title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,ref(v),'r','linewidth',2);
plot(v,dist(v),'g','linewidth',2);
plot(v,u(v),'m','linewidth',2);
legend('Output','target','disturbance','input');
subplot(312);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(313);plot(v,c(:,v)','linewidth',2);title('c_k=0 for OMPC','fontsize',18)

