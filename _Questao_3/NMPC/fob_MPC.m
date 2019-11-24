function J = fob_MPC(xk,mv,uk_1,Hp,Hc,sysi,ysp,q,r,Qa)
%function J = fob_MPC(xk,mv,uk_1,Hp,Hc,sysi,ysp,q,r,Qa)
nu = size(sysi.B,2) ;

% Expansão do vetor u dentro do horizonte de controle/predição
du = [mv;zeros(nu*(Hp-Hc),1)];
%Qa=150*0.001/3600;
% Simulação do modelo
xk = xk + [0.405;0.405; 0.2142;0.2411];
y = simulate_MIMO_incremental(xk,du,uk_1,Hp,sysi,Qa);
%y = y - [0.405;0.405];

Q = diag(repmat(q,1,Hp)); R = diag(repmat(r,1,Hp));

ysp = repmat(ysp,Hp,1);

J = (y-ysp)'*Q*(y-ysp) + du'*R*du;
end
