function [y] = simulate_MIMO_incremental(xk,du,uk_1,Hp,sysi,Qa)

nu = size(sysi.B,2) ;
ny = size(sysi.C,1) ;

h0 = xk(1:2);
cont=1;
    for k = 1:Hp
        
    uk_1 = uk_1 + du(cont:cont+1,1) ;
        
    tspan = ([0 0.1] + (k-1)*[0.1 0.1])';

    [t,ys] = ode45(@(t,h)twotanksODE(t,h,uk_1,Qa),tspan, h0) ;
    
    y(ny*(k-1)+1:ny*k,1) = ys(end,:);

    h0 = ys(end,:);    
        cont=cont+2;
    end

end

