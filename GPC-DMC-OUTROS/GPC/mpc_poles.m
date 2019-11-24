%%%%%  Find GPC closed-loop pole polynomial given control law 
%%%%%  Pr r - Nk y = [1,Dk] Delta u     G=z^{-1}b/a
%%%%%    [No T-filter]
%%%%%   Pc = [I+z^{-1}Dk(z)]*a(z)*Delta(z) + Nk(z)*z^{-1}b(z)
%%%%%
%%%%%  Pc = mpc_poles(a,b,Nk,Dk)
%%%%%
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)
    
function Pc = mpc_poles(a,b,Nk,Dk)

%%% Closed-loop poles without T-filter are given from
%%% Pc = [I+z^{-1}Dk(z)]*a(z)*Delta(z) + Nk(z)*z^{-1}b(z)
    Pc1 = conv([1,Dk],conv(a,[1,-1])); n1 = length(Pc1);
    Pc2 = conv([0,b],Nk); n2 = length(Pc2);
    if n1>n2; Pc2(n1)=0;elseif n2>n1;Pc1(n2)=0;end
    Pc=Pc1+Pc2; 