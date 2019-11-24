function [du, duk] = MPC_twotanks(xmk, Am, Bm, Cm, Hp, Hc, ysp, q, r, du0,uk_1)

sysi.A = Am;
sysi.B = Bm;
sysi.C = Cm;
uk_1 = uk_1+[0.2142;0.2411];
nu = size(sysi.B,2) ;
dumax = ([0.1;0.1;0.1;0.1;0.1;0.1]); dumin = -dumax;
umax = [0.95;0.95];  umin = [0.01; 0.01];
Mtil=[];
Itil=[];
auxM=zeros(nu,Hc*nu);
for in=1:Hc
    auxM=[eye(nu) auxM(:,1:(Hc-1)*nu)];
    Mtil=[Mtil;auxM];
    Itil=[Itil;eye(nu)];
end
%Ain = [Mtil;-Mtil];

du = fmincon(@(du) fob_twotanks_MPC(xmk,du,Hp,Hc,sysi,ysp,q,r),du0,[],[],[],[],dumin,dumax);

%%% COLOCAR RESTRIÇÃO EM U FUNCIONA ATÉ O DISTURBIO EM Q, DEPOIS DÁ
%%% PROBLEMA. NÃO CONSEGUIMOS RESOLVER AINDA.
% du = fmincon(@(du) fob_twotanks_MPC(xmk,du,Hp,Hc,sysi,ysp,q,r),du0,Ain,calculo_Bin(uk_1,umax,umin,Hc,Itil),[],[],dumin,dumax);


duk = du(1:nu) ;
du = [du(nu+1:end);zeros(nu,1)];
end