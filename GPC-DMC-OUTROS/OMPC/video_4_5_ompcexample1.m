%%%%%  Typical data entry required for a SOMPC
%%%%%  SISO state space  example
%%%%%    x(k+1) = A x(k) + B u(k)
%%%%%    y(k) = C x(k) + D u(k)      Note: Assumes D=0
%%%%%
%%%%%   Assumes J = sum x(k+i) Q x(k+1) + u(k+i-1) R u(k+i-1)
%%%%%
%%%%%   and uses  TERMINAL MODE  U = -K x BASED ON Q2,R2
%%%%%
%%%%%       [NOMINAL CASE ONLY]
%%%%%
%%%%%  Illustrates closed-loop simulations with constraint handling
%%%%%  
%%%%%
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%% Model
A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
D=zeros(2,2);

%% Tuning for PERFORMANCE INDEX
Q=C'*C;
R=eye(2);
nc=5;   %%% no. of d.o.f.
%%% parameters for TERMINAL MODE
Q2=Q; R2=R;


