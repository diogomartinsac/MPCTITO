%%%%%   THIS IS A SCRIPT FILE. 
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%
%% demonstrates use of MATLAB code to find GPC law and pole polynomial
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%%% SISO Model and GPC parameters
a=[1 -.9 0.6];
b=[1 0.4];
sizey=1;  %%% siso
ny=20;
nu=2;
lambda=1000;
ny2=60;

%%%% plotting OL vs CL
%%% closed-loop
%%% Set point, disturbance and noise
ref = [zeros(1,25),ones(1,ny2)];
dist=[zeros(1,ny2+25)];
noise = [zeros(1,ny2+25)];
y=cell(1);u=y;Du=y;
for na=1:15;
[y{na},u{na},Du{na},r] = mpc_simulate_noconstraints_withff(b,a,nu,ny,lambda,1,ref,dist,noise,na);
J(na)=norm(r-y{na},2)^2+lambda*norm(Du{na},2)^2;
end

%[y2,u2,Du2,r] = mpc_simulate_noconstraints(b,a,nu,ny,lambda,1,ref,dist,noise);

figure(1); clf reset
v=0:ny;v2=0:ny2;
aa=plot([0 20 20 ny2],[0 0 1 1],'g--'); hold on
for na=1:6;
    %plot(0:ny2-5,u{na}(5:ny2))
    plot(0:ny2-5,y{na}(5:ny2),'color',[na/6,0,1-na/6]);
end
xlim([15,30])
set(aa,'linewidth',2);
c=legend('r','na=1','na=2','na=3','na=4','na=5','na=6');
set(c,'fontsize',14)
title(['n_y =',num2str(ny),' n_u =',num2str(nu),', \lambda = ',num2str(lambda)],'fontsize',18)

figure(2); clf reset
plot(1:15,J);
title('Runtime cost vs  na','fontsize',18)
xlabel('n_a','fontsize',18);
ylabel('J','fontsize',18);
