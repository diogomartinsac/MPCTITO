A =[0.9      0;
    0.2    0.1];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
Q=C'*C;Q2=Q;
R=eye(1); R2=R*10;
nx=2;
nu=1;
K0=dlqr(A,B,Q2,R2);
%%% Horizon 1
nc=4;runtime=30;x0=[3;2];
[J,x,y,u,c,KSOMPC] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime);
figure(1); clf reset
v=2:length(y);
subplot(311);plot(v,y(v),'b','linewidth',2);title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
subplot(312);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(313);
c=[c;zeros(10,size(c,2))];
for k=1:nc+10;
    alp=(nc+10-k)/(nc+10);
    plot(k+1:11+nc,c(1:10+nc-k+1,k)','linewidth',2,'color',[alp,0,1-alp]);hold on
    title('Optimised cfut for SOMPC','fontsize',18)
end
legend('c(k)','c(k+1)','c(k+2)','c(k+4)','c(k+5)')
%%% Investigating feedback dependence on nc
runtime=2;Ksave=[K0];
for nc=1:10;
  [J,x,y,u,c,KSOMPC] = chap4_ompc_simulate(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime);
  Ksave(nc+1,:)=KSOMPC;
end
figure(2); clf reset
plot(0:10,Ksave); title('Parameters of KSOMPC vs n_c','fontsize',18);
xlabel('n_c','fontsize',18);
