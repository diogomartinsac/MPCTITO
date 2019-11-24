function [du, duk] = MPC_twotanks(xmk, Am, Bm, Cm, Hp, Hc, ysp, q, r, du0)

sysi.A = Am;
sysi.B = Bm;
sysi.C = Cm;
nu = size(sysi.B,2) ;

du = fminunc(@(du) fob_twotanks_MPC(xmk,du,Hp,Hc,sysi,ysp,q,r),du0);

duk = du(1:nu) ;
du = [du(nu+1:end);zeros(nu,1)];
end