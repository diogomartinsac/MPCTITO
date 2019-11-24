A =[0.8      0;
    0.2    0.2];
B =[0.2;0.8];
C =[1.9 -1];
Ap=A;Bp=B;Cp=C;
Ap(2,1)=0.3;Ap(2,2)=0.3;Bp(1,1)=0.4;Cp(1,1)=2;
D=0;
Q=C'*C;Q2=Q;
R=1*eye(1); R2=R;
nx=2;
nu=1;
ref=[zeros(1,20),ones(1,50)];
dist=[zeros(1,30),-0.5*ones(1,30)];
noise=dist*0;
nc=2;
na=5;
%%% Horizon 1
runtime=49;x0=[0;0];
[x,y,u,c,z] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);
na=1;
[x2,y2,u2,c2,z2] = chap6_ompc_simulate_withtracking(A,B,C,D,Ap,Bp,Cp,nc,na,Q,R,Q2,R2,x0,runtime,ref,dist,noise);

%%%%% plotting
figure(1); clf reset
subplot(211);
v=2:size(y,2);
plot(v,y(:,v)','b','linewidth',2);title(['OMPC output for {n_c,n_a}=',num2str([nc,5])],'fontsize',18)
hold on; plot(v,ref(:,v)','r','linewidth',2);
%plot(v,dist(:,v)','g','linewidth',2);
plot(v,u(:,v)','m','linewidth',2);
%plot(v,z(nx+1:end,v)','c','linewidth',2)
legend('Output 1','target 1','input 1');
ylim([0 1.6])
subplot(212);
plot(v,y2(:,v)','b','linewidth',2);title(['OMPC output for {n_c,n_a}=',num2str([nc,1])],'fontsize',18)
hold on; plot(v,ref(:,v)','r','linewidth',2);
%plot(v,dist(:,v)','g','linewidth',2);
plot(v,u2(:,v)','m','linewidth',2);
%plot(v,z(nx+1:end,v)','c','linewidth',2)
ylim([0 1.6])


