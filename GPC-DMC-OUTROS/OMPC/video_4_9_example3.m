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
  %R2=R; %to get OMPC
nx=3;
nu=2;
ref=[zeros(2,5),[ones(1,50);zeros(1,20),0.5*ones(1,30)]];
dist=[zeros(2,40),ones(2,15)*0.4];

%%% Horizon 1
nc=5;runtime=50;x0=[0;0;0];
[J,x,y,u,c,KSOMPC] = chap4_ompc_simulateb(A,B,C,D,nc,Q,R,Q2,R2,x0,runtime,ref,dist);

%%%%% plotting
figure(1); clf reset
v=2:size(y,2);
subplot(311);plot(v,y(:,v)','b','linewidth',2);title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,ref(:,v)','r','linewidth',2);
plot(v,dist(:,v)','g','linewidth',2);
plot(v,u(:,v)','m','linewidth',2);
legend('Output 1','Output 2','target 1','target 2','disturbance 1','disturbance 2');

subplot(312);plot(v,J(v),'b','linewidth',2);title('cost is monotonic','fontsize',18)
subplot(313);
c=[c;zeros(10*nu,size(c,2))];
for k=5:nc+10;
    alp=(nc+10-k)/(nc+5);
    for kk=1:nu;
        subplot(3,nu,2*nu+kk);
    plot(k:k+nc+9,c(kk:nu:end,k)','linewidth',2,'color',[alp,0,1-alp]);hold on
    title(['cfut for loop ',num2str(kk)],'fontsize',18)
    xlim([5,20]);
    end
    
end
legend('c(k)','c(k+1)','c(k+2)','c(k+3)','c(k+4)')
