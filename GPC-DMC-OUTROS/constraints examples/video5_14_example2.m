% 
Phi=[0.7,0.4;-0.8,0.9];
A_y = [eye(2,2); -eye(2,2)];
b_y = [10; 10; 10; 10];

[A_S, b_S]=construct_mas(Phi,A_y,b_y,2);

%%% State evolutions
A=Phi;
AA=[eye(2);A;A^2;A^3;A^4;A^5;A^6;A^7;A^8];
x1=AA*[10;10];
x2=AA*[-10;10];
x3=AA*[5;10];
x4=AA*[-4;-10];
x5=AA*[-7;2];
x6=AA*[5;-10];
v=1:2:16;v2=2:2:16;

figure(3); clf reset
Q=Polyhedron(A_y,b_y);
P=Polyhedron(A_S,b_S);  %%% from mpt3 toolbox
plot(Q,'color','y');hold on
plot(P,'color','w'); hold on
legend('Sample constraints','MCAS')
plot(x1(v),x1(v2),'r+-');
plot(x6(v),x6(v2),'r+-');
plot(x2(v),x2(v2),'r+-');
plot(x3(v),x3(v2),'b+-');
plot(x4(v),x4(v2),'b+-');
plot(x5(v),x5(v2),'g+-');
title('MCAS and state evolutions')
axis([-14,14,-14,14])
