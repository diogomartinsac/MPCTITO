%%%% Constructs all  data needed for a SOMPC simulation
%%% into a structure array
%% 
%%%     x(k+1) = A x(k) + B u(k)            x0 is the initial condition
%%%     y(k) = C x(k) + D u(k) + dist       Note: Assumes D=0, dist unknown
%%%
%%  input constraint    umax, umin
%%  state constraints  | Kxmax x | < xmax
%%  Horizon for constraint handling   npred
%%  reference trajectory     ref, r
%%  perturbation to control  c
%%  Number of d.o.f.         nc
%%  Weighting matrices in J  Q, R
%%  Weighting matrices for terminal mode   Q2, R2
%%  measurement noise        noise

%%  simdata=sompc_data(A,B,C,D,Q,R,Q2,R2,nc,umin,umax,Kxmax,xmax,npred)

%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

function simdata=sompc_data(A,B,C,D,Q,R,Q2,R2,nc,umin,umax,Kxmax,xmax,npred)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       Control law is   u = -Knew z + Pr r + c
%%%%
%%%%%      Observor is   z = Ao*z +Bo*u + L*(y + noise - Co*z );        z=[xhat;dhat]
%%%%%
%%%%%      K the underlying control law:  u-uss = -K(x-xss) 
[K,L,Ao,Bo,Co,Do,Knew,Pr,Mx,Mu] = sompc_observor(A,B,C,D,Q,R);

%%% Find predictions based on given control law K
%%%%%   x =  Pc1*c + Pz1*z + Pr1*r + Pd1*d
%%%%%   u =  Pc2*c + Pz2*z + Pr2*r + Pd2*d
[Pc1,Pc2,Pc3,Pz1,Pz2,Pz3,Pr1,Pr2,Pr3] = sompc_predclp(A,B,C,D,Knew,Pr,nc,npred);

%%%%%   Constraints are summarised as 
%%%%%   CC c - dfixed - dx0*[z;r;d] <= 0     d is a known disturbance
%%%%%   Predictions are 
%%%%%   x =  Pc1*c + Pz1*z + Pr1*r + Pd1*d
%%%%%   u =  Pc2*c + Pz2*z + Pr2*r + Pd2*d
[CC,dfixed,dx0] = sompc_constraints(Pc1,Pc2,Pz1,Pz2,Pr1,Pr2,umin,umax,Kxmax,xmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Performance is based on Q2, R2 with predictions based on Q,R:
%%%  measured as
%%  x SX x + x SXC c + c SC c
%[SX,SC,SXC]=sompc_cost(A,B,C,D,Q,R,Q2,R2,nc);
[SX,SC,SXC]=sompc_cost2(A,B,C,D,Q,R,Q2,R2,nc);


%keyboard
simdata=struct('L',L,'K',K,'Knew',Knew,'A',A,'B',B,'C',C,'D',D,'nc',nc,'CC',CC,'dfixed',dfixed,...
    'dx0',dx0,'SX',SX,'SXC',SXC,'SC',SC,'Ao',Ao,'Bo',Bo,'Co',Co,'Do',Do,'Pr',Pr,'Mx',Mx,'Mu',Mu,...
    'Q2',Q2,'R2',R2);


