function J = fob_MPC(xk,mv,Hp,Hc,sysi,ysp,q,r)

nu = size(sysi.B,2) ;

% Expans�o do vetor u dentro do horizonte de controle/predi��o
du = [mv;zeros(nu*(Hp-Hc),1)];

% Simula��o do modelo
y = simulate_MPC(xk,du,Hp,sysi);

Q = diag(repmat(q,1,Hp)); R = diag(repmat(r,1,Hp));

ysp = repmat(ysp,Hp,1);

J = (y-ysp)'*Q*(y-ysp) + du'*R*du;
end
