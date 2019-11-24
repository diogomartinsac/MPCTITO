% construct double integrator with uncertainty
Phi=[0.7,0.4,0.2;-0.8,0.9,0.3;-0.5,0,0.75]*0.9;
Phi(4,4)=1;
A_y2 = [[eye(3,3); -eye(3,3)],zeros(6,1)];
b_y2 = [1; 1;1; 2; 3;4];
A_y1=[zeros(2,3),[1;-1]];
b_y1=[1;1];

[A_S, b_S]=construct_mas_tracking(Phi,A_y1,b_y1,A_y2,b_y2,2);
