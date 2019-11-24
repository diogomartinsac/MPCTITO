%%%  Computes and simulates a GPC control law using an independent model
%%%  approach to prediction
%%%  Assumes no advance knowledge of target is used
%%%
%%%%%  [y,u,Du,r] = immpc_tf_simulate_nonconstraints(B,A,nu,ny,Wu,Wy,ref,dist,noise)
%              y, u, Du, r are dimensionally compatible 
% closed-loop outputs/inputs/input increments and supplied set-point and disturbance
%
% MFD model     Ay(k) = Bu(k-1) + dist
%
% ny is output horizon
% nu is the input horizon
% Wu is the diagonal control weighting 
% Wy is the diagonal output weighting
% sizey no. outputs and inputs (assumed square)
% dist,noise are the disturbance and noise signals
% ref is the reference signal
% Dumax is a vector of limits on input increments (assumed symetric)
% umax, umin are vectors of limits on the inputs
% ymax, ymin are vectors of limits on outputs
%
% [y,u,Du,r,d] = imgpc_tf_simulate_constraints(B,A,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise)
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

function [y,u,Du,r] = imgpc_tf_simulate_constraints(B,A,nu,ny,Wu,Wy,Dumax,umax,umin,ymax,ymin,ref,dist,noise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error checks
sizey = size(A,1);
if size(B,2)==sizey;B=[B,zeros(sizey,sizey)];end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Find prediction matrices 
%%%%    yfut = H *Dufut + P*Dupast + Q*ypast+L*offset
[H,P,Q,L] = imgpc_tf_predmat(A,B,ny);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Define constraint matrices
%%%%%%   CC*Du(future) - dd - du*u(k-1))-ddu*Dupast-dy*ypast-ddd*dk <= 0
[CC,dd,du,ddu,dy,ddd]  = imgpc_constraints_tfmodel(Dumax,umax,umin,ymax,ymin,sizey,nu,H(:,1:sizey*nu),P,Q,L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Find control law and parameters of the cost function
%%%%   Dufut = Pr*rfut - Dk*Dupast - Nk*ypast - Mk*offset
%%%%    J = Dufut'*S*Dufut + Dufut'*2X*[Dupast;ypast;rfut;offset]
[Nk,Dk,Pr,Mk,S,X] = imgpc_tf_law(H,P,Q,L,nu,Wu,Wy,sizey);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Set up simulation parameters
nNk = size(Nk,2)/sizey;
nDk = size(Dk,2)/sizey;
init = max([nNk,nDk])+2;
y = zeros(sizey,init);
ym=y;
u = y;
Du = u;
r = u;
d=u;
offset=d;
runtime = size(ref,2);
opt = optimset('quadprog');
opt.Diagnostics='off';    %%%%% Switches of unwanted MATLAB displays
opt.LargeScale='off';     %%%%% However no warning of infeasibility
opt.Display='off';
opt.Algorithm='active-set';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf reset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Closed-loop simulation
for i=init:runtime-2;

%%% Update unconstrained control law
d(1:sizey,i+1)=dist(:,i+1);
ympast = ym(:, i:-1:i+1-nNk);
offset(:,i)=y(:,i)-ym(:,i)+noise(:,i); %% noise does not change actual output
Dupast = Du(:, i-1:-1:i-nDk) ;
upast = u(:, i-1);
rfut = ref(:,i+1); 

%%%%%%% Unconstrained law 
Dufast(:,i) = Pr*rfut - Nk*ympast(:) - Dk*Dupast(:) - Mk*offset(:,i);
ufast(:,i)=u(:,i-1)+Dufast(1:sizey,i);

% Form constraint matrices and solve constrained optimisation
%%%%%%   CC*Du(future) - dd - du*u(k-1))-ddu*Dupast-dy*ypast <= 0
%  J = Dufut'S*Dufut+2*X*[Dupast(:);ypast(:);rfut(:);offset(:)]*Dufut
dt = dd+du*upast(:)+ddu*Dupast(:)+dy*ympast(:)+ddd*offset(:,i);
f=X*[Dupast(:);ympast(:);rfut(:);offset(:,i)];
[Duqp,fval,exitflag] = quadprog(S,f,CC,dt,[],[],[],[],[],opt);
if exitflag==-2; %%% output constraints causing a problem so ignore
    disp(['Ignoring output constraints at sample ',num2str(i)]);
  vec=1:4*nu*sizey;
  [Duqp,fval,exitflag] = quadprog(S,f,CC(vec,:),dt(vec),[],[],[],[],[],opt);
end

Du(:,i) = Duqp(1:sizey);
u(:,i) = u(:,i-1)+Du(:,i);


%%% Simulate the process
upast2 = u(:,i:-1:i-nDk);
ypast2 = y(:, i:-1:i+2-nNk);
y(:,i+1) = -A(:,sizey+1:nNk*sizey)*ypast2(:) + B*[upast2(:)] + d(:,i+1);
r(:,i+1) = ref(:,i+1);

%%% Simulate the model (for now same parameters but no disturbance/noise)
umpast2 = u(:,i:-1:i-nDk);
ympast2 = ym(:, i:-1:i+2-nNk);
ym(:,i+1) = -A(:,sizey+1:nNk*sizey)*ympast2(:) + B*[umpast2(:)];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Produce a neat plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Ensure all outputs are dimensionally compatible
u(:,i+1) = u(:,i);
Du(:,i+1) = Du(:,i)*0;
noise = noise(:,1:i+1);
ufast(:,i+1)=ufast(:,i);
Dufast(:,i+1)=Dufast(:,i)*0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


time=0:size(u,2)-1;
for i=1:sizey;
    figure(i);clf reset
    plotall(y(i,:),r(i,:),u(i,:),ufast(i,:),Du(i,:),Dufast(i,:),d(i,:),noise(i,:),umax(i),umin(i),Dumax(i),ymax(i),ymin(i),time,i);
end

disp('*******************************************************************************');
disp(['***    For GPC there are ',num2str(sizey),' figures beginning at figure 1   ***']);
disp('*******************************************************************************');


%%%%% Function to do plotting in the MIMO case and 

function plotall(y,r,u,ufast,Du,Dufast,d,noise,umax,umin,Dumax,ymax,ymin,time,loop)

uupper = [umax,umax]';
ulower = [umin,umin]';
yupper = [ymax,ymax]';
ylower = [ymin,ymin]';
Dulim = [Dumax,Dumax]';
time2 = [0,time(end)];
rangeu = (max(umax)-min(umin))/10;
rangey = (max(ymax)-min(ymin))/10;
ranged = (max(max([d,noise]))-min(min([d,noise])))/20;if ranged==0;ranged=1;end

subplot(221);plot(time,y','-',time,r','--',time2,yupper,'--',time2,ylower,'--');
axis([time2,ymin-rangey,ymax+rangey]);
xlabel(['IMGPC - Outputs and set-point in loop ',num2str(loop)]);
subplot(222);plot(time,Du','-',time, Dufast,':',time2,Dulim,'--',time2,-Dulim,'--');
axis([time2,min(-Dumax)-rangeu,max(Dumax)+rangeu]);
legend('Constrained','Unconstrained');
xlabel(['IMGPC - Input increments in loop ',num2str(loop)]);
subplot(223);plot(time,u','-',time,ufast,':',time2,uupper,'--',time2,ulower,'--');
axis([time2,min(umin)-rangeu,max(umax)+rangeu]);
xlabel(['IMGPC - Inputs in loop ',num2str(loop)]);
subplot(224);plot(time,d','b',time,noise,'g');
axis([time2,min(min([d,noise]))-ranged,max(max([d,noise]))+ranged]);
xlabel(['IMGPC - Disturbance/noise in loop ',num2str(loop)]);
