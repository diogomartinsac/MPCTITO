A =[0.8      0.4 1;
    -0.2    0.7 -0.2;1.2 0.3 0.2];
B =[0.4;0.8;0];
C =[1.9 2 -0.2];
D=0;
Q=C'*C;Q2=Q;
R=2*eye(1); R2=R*10;
nx=3;
nu=1;
K0=dlqr(A,B,Q2,R2);
ref=[zeros(1,10),ones(1,50)];
dist=[zeros(1,30),-0.5*ones(1,30)];

%%% Horizon 1
nc=2;runtime=49;x0=[0;0;0];
[J,x,y,u,c,KSOMPC] = chap4_ompc_simulateb(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist);
figure(1); clf reset
v=2:length(y);
subplot(311);plot(v,y(v),'b','linewidth',2);title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,ref(v),'r','linewidth',2);
plot(v,dist(v),'g','linewidth',2);
plot(v,u(v),'m','linewidth',2);
legend('Output','target','disturbance','input');
subplot(312);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(313);
c=[c;zeros(20,size(c,2))];
for k=9:nc+12;
    alp=(nc+12-k)/(3+nc);
    plot(k+1:21+nc,c(1:20+nc-k+1,k)','linewidth',2,'color',[alp,0,1-alp]);hold on
    title('Optimised cfut for SOMPC','fontsize',18)
end
legend('c(k)','c(k+1)','c(k+2)','c(k+4)','c(k+5)')

