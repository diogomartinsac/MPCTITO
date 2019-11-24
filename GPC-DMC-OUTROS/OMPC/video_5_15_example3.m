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
nc=3;

%%%%% constraints
umin=-1;     %%% umin<u<umax
umax=1.35;
Kxmax=[1 0.2;-0.1 0.4;-1,-0.2;0.1,-0.4];
xmax=[4;4;0.8;2.5];
rdmax=20;
rdmin=-20;


%%%  %%%% Check that the target is OK
%%%%  means all gaps should be negative
gaps = feasibletargetcheck(A,B,C,umin,umax, Kxmax, xmax,rdmax,rdmin);
if any(gaps>0); disp('infeasible limits on targets');end
    
%%% 
ref=[zeros(1,10),ones(1,85)*3];
dist=[zeros(1,30),-ones(1,65)*0.4];

%%% initial condition
x0=[1.3;-2]*0;

%%% Horizon 1
nc=2;runtime=29;

[J,x,y,u,c,Ksompc,F,t,G1,f1] = ompc_simulate_efficient(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist,umin,umax,Kxmax,xmax,rdmax,rdmin);

figure(1); clf reset
v=2:length(y);
subplot(221);plot(v,y(v),'b','linewidth',2);hold on
plot(v,u(v),'m','linewidth',2);
legend('output','constrained u')
plot([0,runtime],[umax,umax],'m--',[0,runtime],[umin,umin],'m--');
title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(222);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(223);plot(v,c(:,v)','linewidth',2);title('c_k for OMPC','fontsize',18)
subplot(224);plot(x(1,:),x(2,:),'r','linewidth',2);title('State plane','fontsize',18);



%%  Mx +Nc+V*(r-d) <= t 
P=Polyhedron(F,t);
R=Polyhedron(G1,f1);
Qp = projection(P, [5]); %%% from mpt3 toolbox
figure(4); clf reset
plot(Qp,'color','y'); hold on
title('Achieveable values of r-d')


