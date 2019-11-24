%%%%%  Typical data entry required for a 
%%%%%  SISO transfer function example
%%%%%    A(z) y  =  z^{-1}B(z) u
%%%%%
%%%%%  Illustrates closed-loop GPC simulations with NO constraint handling
%%%%%
%%% MIMO Model and GPC parameters
A=[];B=[];
A(1,1:3:12) = poly([.5,.8,-.2]);
A(2,2:3:12) = poly([-.2,.9,.5]);
A(3,3:3:12) = poly([-.1,.4,.6]);
A(1,5:3:12) = [ -.2 .1 .02];
A(2,4:3:12) = [ .4 0 -.1];
A(3,4:3:12) = [0 .2 .2];
A(2,6:3:12) = [-1 .1 .3];
A(:,4:12) = A(:,4:12)*.8;

B = [.5 0.2 -.5 1 2 1;2 0 .3 -.8 .6 .5;0 .9 -.4 1 .3 .5];
sizey=3;
Tfilt=[eye(3),-0.8*eye(3)];

%% Tuning parameters
Wu =eye(3); % matrix of input weights
Wy=eye(3);  % matrix of output weights
ny=9;  % prediction horizon
nu=4;   % input horizon
Dumax=[0.1;0.1;0.2];
umax=[1.2;0.5;1];
umin=[-0.3;-0.4;-1.9];
ymin=[-0.2;-0.3;-0.2];
ymax=[1.1;0.05;1];

%%% Set point, disturbance and noise
ref = [zeros(1,10),ones(1,30);zeros(1,10),-0.2*ones(1,30)/2;zeros(1,10),ones(1,30)*.7];
dist=[zeros(3,10),ones(3,80)*0.2]*0; 
noise = [zeros(3,25),randn(3,40)*0.01]*0;

%%%% Add a disturbance
dist=[zeros(3,20),ones(3,70)*0.2]; 
[y,u,Du,r] = mpc_simulate_tfiltoutputconstraints(B,A,Tfilt,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);

[y,u,Du,r] = imgpc_tf_simulate_constraints(B,A,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise);
