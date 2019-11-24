%%% Contours for a performance index
%%%  J=0.5x'Sx+b'*x+c
b=[1.5 20]'; 
c=3;
S=[0.8 0.8;0.8 5.8];
M=[-1 0;0 -1;1 0;0 1];
d=[1;1;1;1];
opt.Algorithm='active-set';

xunconstrained=-inv(S)*b;
xopt=quadprog(S,b,M,d,[],[],[],[],[],opt);

x=linspace(-5,5,300);
y=x-3;

for k=1:length(x);
    for kk=1:length(y);
        z=[x(k);y(kk)];
        J(k,kk)=0.5*z'*S*z+z'*b+c;
    end
end
figure(1); clf reset;
contour(x,y,J',12);hold on
plot(xopt(1),xopt(2),'o','markersize',20)
plot(xunconstrained(1),xunconstrained(2),'b+','markersize',20)

grid

plot([1 1],[-8 5],'m--','linewidth',2);
plot([-1 -1],[-8 5],'m--','linewidth',2);
plot([-5,5],[1 ,1],'m--','linewidth',2);
plot([-5,5],[-1,-1],'m--','linewidth',2);

