% construct double integrator with uncertainty
Phi=[0.7,0.4,0.2;-0.8,0.9,0.3;-0.5,0,0.75]*0.9;
A_y = [eye(3,3); -eye(3,3)];
b_y = [1; 1;1; 2; 3;4];

[A_S, b_S]=construct_mas(Phi,A_y,b_y,2);


figure(3); clf reset
Q=Polyhedron(A_y,b_y);
P=Polyhedron(A_S,b_S);  %%% from mpt3 toolbox
plot(P,'color','c'); 
title('MCAS ')
figure(4); clf reset
plot(Q,'color','y'); 
title('Sample constraints')
