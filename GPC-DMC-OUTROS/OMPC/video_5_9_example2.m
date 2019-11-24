A=[0.8,0.1;-0.2,0.9];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
%%% Q, R for J   
%%% Q2, R2 for terminal mode
Q=C'*C;Q2=Q;
R=0.1*eye(1); R2=R;
nx=2;
nu=1;

%%%%% constraints
umin=-1;     %%% umin<u<umax
umax=1;
Kxmax=[1 0.2;-0.1 0.4;-1,-0.2;0.1,-0.4];
xmax=[1;2;0.8;1.5];

%%% initial condition
x0=[0;4.5];
%x0=[-1.5;4.5];
%x0=[1.3;-2];

%%% Horizon 1
nc=4;runtime=14;
[J,x,y,u,c,KSOMPC] = chap5_ompc_simulate_constraints(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,umin,umax,Kxmax,xmax);
%%% unconstrained
[Ju,xu,yu,uu,cu,KSOMPC] = chap5_ompc_simulate_constraints(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,umin*100,umax*100,Kxmax,xmax*100);
figure(1); clf reset
v=2:length(y);
subplot(221);plot(v,y(v),'b','linewidth',2);hold on
plot(v,u(v),'m','linewidth',2);
plot(v,uu(v),'g','linewidth',2);
legend('output','constrained u','unconstrained u')
plot([0,runtime],[umax,umax],'m--',[0,runtime],[umin,umin],'m--');
title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(222);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(223);plot(v,c(:,v)','linewidth',2);title('c_k for OMPC','fontsize',18)
subplot(224);plot(x(1,:),x(2,:),'r','linewidth',2);title('State plane','fontsize',18);

%%%% Data to plot notional constraints on x
xx=[-2,2];
y1=[xmax-Kxmax(:,1)*xx(1)]./Kxmax(:,2);
y2=[xmax-Kxmax(:,1)*xx(2)]./Kxmax(:,2);
figure(2); clf reset
plot(x(1,:),x(2,:),'r+','linewidth',2);hold on
plot(xu(1,:),xu(2,:),'bo','linewidth',2);
legend('Constrained x(k)','unconstrained x(k)');
plot(xx,[y1,y2]','k--'); 
ylim([-5,5])
xlim([-2,1.5])
title('Phase trajectory and constraints at each sample')

%%%% Project a set to a lower dimension
%%  Mx +Nc <= t  ==>  Qx<=E
load settest
P=Polyhedron([M,N],t);  %%% from mpt3 toolbox
P2=Polyhedron(M,t);
Qp = projection(P, [1,2]); %%% from mpt3 toolbox
plot(Qp,'color','w',P2,'color','y'); hold on
ylim([-5,5])
xlim([-2,1.5])
plot(x(1,:),x(2,:),'r+','linewidth',2);hold on
plot(xu(1,:),xu(2,:),'bo','linewidth',2);

