% 
Phi=[0.7,0.4,0;-0.8,0.9,0;0 0 1];

A_y2 = [[eye(2,2); -eye(2,2)],zeros(4,1)];
b_y2 = [10; 10; 10; 10];
A_y1=[zeros(2,2),[1;-1]];
b_y1=[1;1];

[A_S, b_S]=construct_mas_tracking(Phi,A_y1,b_y1,A_y2,b_y2,2);

%%% State evolutions
A=Phi;
AA=[eye(3);A;A^2;A^3;A^4;A^5;A^6;A^7;A^8];
x1=AA*[10;10;0.5];
x2=AA*[-10;10;0.2];
x3=AA*[5;10;0.3];
x4=AA*[-4;-10;-0.5];
x5=AA*[-7;2;-0.6];
x6=AA*[5;-10;1];
v=1:3:27;v2=2:3:27;v3=3:3:27;

figure(3); clf reset
Q=Polyhedron([A_y1;A_y2],[b_y1;b_y2]);
P=Polyhedron(A_S,b_S);  %%% from mpt3 toolbox
%plot(Q,'color','y');hold on
plot(P,'color','w'); hold on
legend('MCAS')
plot3(x1(v),x1(v2),x1(v3),'r+-');
plot3(x6(v),x6(v2),x6(v3),'r+-');
plot3(x2(v),x2(v2),x2(v3),'r+-');
plot3(x3(v),x3(v2),x3(v3),'b+-');
plot3(x4(v),x4(v2),x4(v3),'b+-');
plot3(x5(v),x5(v2),x5(v3),'g+-');
title('MCAS and state evolutions')
axis([-14,14,-14,14])
