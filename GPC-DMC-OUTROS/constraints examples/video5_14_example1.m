% 
Phi=[0.8,0.1;-0.2,0.9];
A_y = [eye(2,2); -eye(2,2)];
b_y = [10; 10; 10; 10];

[A_S, b_S]=construct_mas(Phi,A_y,b_y,2);

figure(3); clf reset
P=Polyhedron(A_S,b_S);  %%% from mpt3 toolbox
plot(P,'color','w'); hold on
title('MCAS')
