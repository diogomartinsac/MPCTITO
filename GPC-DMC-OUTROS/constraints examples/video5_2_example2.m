%%% Contours for a performance index
%%%  J=0.5x'Sx+b'*x+c
b=[1 1]'; 
b=[1 -5]'; %% move centre
c=3;
S=[0.8 0.6;0.6 1.8];
M=[-1 0;0 -1;-1 -1];
d=[-1;-1;-3];
opt.Algorithm='active-set';

xunconstrained=-inv(S)*b
xopt=quadprog(S,b,M,d,[],[],[],[],[],opt)

x=linspace(-4,4,300);
y=x;
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

plot([1 1],[-6,4],'m--','linewidth',2);
plot([-4,4],[1 ,1],'m--','linewidth',2);
plot([-3,4],[6,-1],'m--','linewidth',2);

