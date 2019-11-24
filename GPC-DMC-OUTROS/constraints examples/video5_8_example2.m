%%%%  Compute and plot MAS and state evolutions
%%% Process is x(k+1)=Ax(k)
%%% Constraints are Cx <=f

%%%% 2-state case
A=[0.8,0.1;-0.2,0.9];
C=[1 0.2;-0.1 0.4;-1,-0.2;0.1,-0.4];
f=[1;2;0.8;1.5];
[F,t]=findmas(A,C,f)

%%%% Data to plot set
x=[-2,2];
y1=[f-C(:,1)*x(1)]./C(:,2);
y2=[f-C(:,1)*x(2)]./C(:,2);

%%% Possible initial points based on initial constraints only
x1=inv(C(1:2,:))*f(1:2);
x2=inv(C([1,4],:))*f([1,4]);
x3=inv(C(2:3,:))*f(2:3);
x4=inv(C(3:4,:))*f(3:4);

%%% data to find state evolution
AA=[eye(2);A;A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;A^11;A^12;A^13];
x1=AA*x1;
x2=AA*x2;
x3=AA*x3;
x4=AA*x4;
vec=1:2:16;

%%% data to find state evolution from inside MAS
AA=[eye(2);A;A^2;A^3;A^4;A^5;A^6;A^7;A^8;A^9;A^10;A^11;A^12;A^13];
xx1=AA*[0.6;-2.5];
xx2=AA*[1.3;-1.5];
xx3=AA*[0.6;1.95];

%%%% Shows that state evolution voilates constraints
figure(1); clf reset
plot(x,[y1,y2]','k--');
hold on
plot(x1(vec),x1(vec+1),'r',x2(vec),x2(vec+1),'g',x3(vec),x3(vec+1),'b',x4(vec),x4(vec+1),'m');
text(-1,0,'Just taking initial constraints, simulations infeasible')
ylim([-6,6])

%%%%%% Plot all inequalities in MAS
%%%% Data to plot set
x=[-2,2];
y1=[t-F(:,1)*x(1)]./F(:,2);
y2=[t-F(:,1)*x(2)]./F(:,2);
figure(2); clf reset
plot(x,[y1,y2]','k--');
hold on
plot(x1(vec),x1(vec+1),'c',x2(vec),x2(vec+1),'c',x3(vec),x3(vec+1),'c',x4(vec),x4(vec+1),'c');
plot(xx1(vec),xx1(vec+1),'r',xx2(vec),xx2(vec+1),'g',xx3(vec),xx3(vec+1),'b');
%text(-1,0,'Just taking initial constraints, simulations infeasible')
ylim([-4,4])
xlim([-1.5,1.5])



