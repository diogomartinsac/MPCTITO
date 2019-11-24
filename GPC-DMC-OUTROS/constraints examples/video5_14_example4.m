% 
Phi=[0.8,0.1,0;-0.2,0.9,0;0 0 1];
A_y2 = [[eye(2,2); -eye(2,2)],zeros(4,1)];
b_y2 = [10; 10; 10; 10];
A_y1=[zeros(2,2),[1;-1]];
b_y1=[1;1];

[A_S, b_S]=construct_mas_tracking(Phi,A_y1,b_y1,A_y2,b_y2,2);

figure(3); clf reset
P=Polyhedron(A_S,b_S);  %%% from mpt3 toolbox
plot(P,'color','w'); hold on
title('MCAS')
