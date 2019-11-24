function [y] = simulate_twotanks_MPC(xk,du,Hp,sysi)

nu = size(sysi.B,2) ;
ny = size(sysi.C,1) ;

    for k = 1:Hp
        xk = sysi.A*xk + sysi.B*du(nu*(k-1)+1:nu*k);
        y(ny*(k-1)+1:ny*k,1) = sysi.C*xk;
    end

end

