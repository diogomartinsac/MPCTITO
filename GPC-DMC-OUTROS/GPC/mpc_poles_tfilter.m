%%%%%  Find GPC closed-loop pole polynomial given control law 
%%%%%  Pr r - Nkt y/T = [1,Dkt] Delta u/T     G=z^{-1}b/a
%%%%%    [Assumes a T-filter]
%%%%%
%%% Pc = [Tfilt+z^{-1}Dkt(z)]*A(z)*Delta(z) + Nkt(z)*z^{-1}B(z)
%%%%%
%%%%%  Pc = mpc_poles_tfilter(a,b,Nkt,Dkt,T)
%%%%%
%%  
%% Author: J.A. Rossiter  (email: J.A.Rossiter@shef.ac.uk)
    
function [Pct,TD] = mpc_poles_tfilter(a,b,Nkt,Dkt,Tfilt)

%%% Closed-loop poles without T-filter are given from
%%% Pc = [Tfilt+z^{-1}Dkt(z)]*A(z)*Delta(z) + Nkt(z)*z^{-1}B(z)
    
    TD=Tfilt;if length(Dkt)>length(Tfilt)-1; TD(1,length(Dkt)+1)=0;end
    TD(2:length(Dkt)+1) = TD(2:length(Dkt)+1)+Dkt;
    Pc1 = conv(TD,conv(a,[1,-1])); n1 = length(Pc1);
    Pc2 = conv([0,b],Nkt); n2 = length(Pc2);
    if n1>n2; Pc2(n1)=0;elseif n2>n1;Pc1(n2)=0;end
    Pct=Pc1+Pc2; 

   