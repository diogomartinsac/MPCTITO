A =[0.9      0;
    0.2    0.1];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
Q=C'*C;Q2=Q;
R=eye(1); R2=R;
nx=2;
nu=1;

%%% Horizon 1
nc=2;runtime=30;x0=[3;2];
[J,x,y,u,c,KSOMPC] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime);
figure(1); clf reset
v=2:length(y);
subplot(311);plot(v,y(v),'b','linewidth',2);title(['OMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(312);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(313);plot(v,c(:,v)','linewidth',2);title('c_k=0 for OMPC','fontsize',18)

%%% Investigating feedback dependence on nc
runtime=2;Ksave=[];
for nc=1:10;
  [J,x,y,u,c,KSOMPC] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime);
  Ksave(nc,:)=KSOMPC;
end
figure(2); clf reset
plot(Ksave); title('Parameters of KSOMPC vs n_c','fontsize',18);
xlabel('n_c','fontsize',18);
