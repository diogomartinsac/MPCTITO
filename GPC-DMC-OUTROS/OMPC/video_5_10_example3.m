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

%%%%% constraints
umin=[-1;-2] ;    %%% umin<u<umax
umax=[2;1.5];
Kxmax=[1 0.2 0;-0.1 0.4 0;-1,-0.2 0;0.1,-0.4 0;0 0 1;0 0 -1];
xmax=[4;4;2.5;0.5;2;2];


%%% initial condition
x0=[0;0;0];

%%% Horizon 1
nc=5;runtime=29;

[J,x,y,u,c,Ksompc,F,t,M1] = chap5_ompc_simulate_constraintsb(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,umin,umax,Kxmax,xmax);

figure(1); clf reset
v=2:length(y);
subplot(221);plot(v,y(:,v)','b','linewidth',2);hold on
plot(v,u(:,v)','m','linewidth',2);
plot([0,runtime],[umax(1),umax(1)],'m--',[0,runtime],[umin(1),umin(1)],'m--');
plot([0,runtime],[umax(2),umax(2)],'c--',[0,runtime],[umin(2),umin(2)],'c--');

title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(222);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(223);plot(v,c(:,v)','linewidth',2);title('c_k for SOMPC','fontsize',18)


