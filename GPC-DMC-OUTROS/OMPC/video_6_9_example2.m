A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
Ap=A;Bp=B;Cp=C;
D=zeros(2,2);
Q=C'*C;Q2=Q;
R=eye(2);R2=R;
nx=3;
nu=2;
na=5;
nc=5;
ref=[zeros(2,8),[ones(1,100);zeros(1,10),0.5*ones(1,90)]];
dist=[zeros(2,70),ones(2,35)*0.0];
noise=dist*0;

%%% Horizon 1
runtime=29;x0=[0;0;0];
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
legend('Output 1','Output 2','target 1','target 2','input 1','input 2');
subplot(212)
plot(v,y2(:,v)','b','linewidth',2);title(['OMPC output for {n_c,n_a}=',num2str([nc,1])],'fontsize',18)
hold on; plot(v,ref(:,v)','r','linewidth',2);
%plot(v,dist(:,v)','g','linewidth',2);
plot(v,u2(:,v)','m','linewidth',2);
%plot(v,z(nx+1:end,v)','c','linewidth',2)


