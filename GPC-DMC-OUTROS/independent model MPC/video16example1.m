%%%%%  Typical data entry required for a 
%%%%%  SISO state space  example
%%%%%    x(k+1) = A x(k) + B u(k)
%%%%%    y(k) = C x(k) + D u(k)     Note: Assumes D=0
%%%%%
%%%%%  Illustrates open-loop prediction
%%%%   yfut = P*x + H*ufut + L*offset   
%%%%   [offset = y(process) - y(model)]
%%%%%  
%%%%%   THIS IS A SCRIPT FILE. CREATES ITS OWN DATA AS REQUIRED
%%%%%   EDIT THIS FILE TO ENTER YOUR OWN MODELS, ETC.
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)

%% SISO Model
A=0.9;
B=1;
C=1;
D=0;
ny=5;

[H,P,L] = imgpc_predmat(A,B,C,D,ny)
[Hx,Px] = imgpc_predmat(A,B,1,D,ny)


%%%% MIMO model
A =[0.9146         0    0.0405;
    0.1665    0.1353    0.0058;
         0         0    0.1353];
B =[0.0544   -0.0757;
    0.0053    0.1477;
    0.8647         0];
C =[1.7993   13.2160         0;
    0.8233         0         0];
D=zeros(2,2);
ny=4;

[H,P,L] = imgpc_predmat(A,B,C,D,ny)
[Hx,Px] = imgpc_predmat(A,B,eye(3),D,ny)


