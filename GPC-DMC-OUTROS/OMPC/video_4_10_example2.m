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
R=eye(2);R2=0.1*R;
  %R2=R; %to get OMPC
nx=3;
nu=2;
ref=[zeros(2,5),[ones(1,60);zeros(1,20),0.5*ones(1,40)]];
dist=[zeros(2,40),ones(2,25)*0.4];
noise=dist*0;
%%% Horizon 1
nc=5;runtime=59;x0=[0;0;0];
[J,x,y,u,c,KSOMPC,z,ym] = chap4_ompc_simulatec(A,B,C,D,Ap,Bp,Cp,nc,Q,R,Q2,R2,x0,runtime,ref,dist,noise);

%%%%% plotting
figure(1); clf reset
v=2:size(y,2);
plot(v,y(:,v)','b','linewidth',2);title(['SOMPC output for n_c=',num2str(nc)],'fontsize',18)
hold on; plot(v,ref(:,v)','r','linewidth',2);
plot(v,dist(:,v)','g','linewidth',2);
plot(v,u(:,v)','m','linewidth',2);
plot(v,ym(:,v)','k--','linewidth',2);
legend('Output 1','Output 2','target 1','target 2','disturbance 1',...
    'disturbance 2','input 1','input 2','model');



