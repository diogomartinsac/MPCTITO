A =[0.9      0;
    0.2    0.1];
B =[0.1;0.8]*1;
C =[1.9 -1];
D=0;
%%% Q, R for J   
%%% Q2, R2 for terminal mode
Q=C'*C;Q2=Q;
R=0.1*eye(1); R2=R;
nx=2;
nu=1;

%%%%% constraints
umin=-0.3;     %%% umin<u<umax
umax=0.2;
Kxmax=[eye(2);-eye(2)];  %%% Kxmax*x < xmax
xmax=[3;3;2;2];

%%% Horizon 1
nc=2;runtime=29;x0=[3;3];
[J,x,y,u,c,KSOMPC] = chap5_ompc_simulate_constraints(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,umin,umax,Kxmax,xmax);
figure(1); clf reset
v=2:length(y);
subplot(221);
%plot(v,y(v),'b','linewidth',2);hold on
plot(v,u(v),'m','linewidth',2);hold on
plot([0,runtime],[umax,umax],'m--',[0,runtime],[umin,umin],'m--');
title(['OMPC input for n_c=',num2str(nc)],'fontsize',18)
subplot(222);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(223);plot(v,c(:,v)','linewidth',2);title('c_k=0 for OMPC','fontsize',18)
subplot(224);plot(x(1,:),x(2,:),'r','linewidth',2);title('State plane','fontsize',18);
