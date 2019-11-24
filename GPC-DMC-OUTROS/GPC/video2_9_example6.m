%%%%%  Typical data entry required for a 
%%%%%  SISO transfer function example
%%%%%    A(z) y  =  z^{-1}B(z) u
%%%%%
%%%%%  Illustrates closed-loop GPC sensitivity for 
%%%%%    different choices of T-filter
%%%%%
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%% Model
a=[1 -1.2 0.32]; 
b=[1,.3];
Tfilt1=[1 -1.6 0.64];
Tfilt2=[1 -0.8];
sizey=1;
Wu =1; % input weights
Wy=1;  % output weights
ny=15;  % prediction horizon
nu=3;   % input horizon

%%%% calculate sensitivity (as relative gains only)
[mag1,mag2,w,mag1y,mag2y]=mpc_sensitivity(a,b,Tfilt1,Wu,Wy,ny,nu);
[mag3,mag4,w2,mag3y,mag4y]=mpc_sensitivity(a,b,Tfilt2,Wu,Wy,ny,nu);

%%%%% Plot input sensitivity
figure(1); clf reset;
loglog(w,mag1(:),'b-',w,mag2(:),'r:');hold on
loglog(w2,mag4(:),'g:')
a=legend('1',num2str(Tfilt1),num2str(Tfilt2));set(a,'fontsize',20);
title('Input sensitivity to measurement noise for different T','fontsize',20)

%%%%% Plot output sensitivity
figure(2); clf reset;
loglog(w,mag1y(:),'b-',w,mag2y(:),'r:');hold on
loglog(w2,mag4y(:),'g:')
a=legend('1',num2str(Tfilt1),num2str(Tfilt2));set(a,'fontsize',20);
title('Output sensitivity to disturbances for different T','fontsize',20)


