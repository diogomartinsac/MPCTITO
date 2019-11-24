%%% Contours for a performance index
%%%  J=x'Sx+x'b+c
b=[4 6]'; 
b=[5 0]'; 
%b=[1 3]';
c=3;
S=[1.2 -0.2;-0.2 0.8];
x=linspace(-4,4,300);
y=x-2;
xopt=-inv(S)*b/2;

for k=1:length(x);
    for kk=1:length(y);
        z=[x(k);y(kk)];
        J(k,kk)=z'*S*z+z'*b+c;
        J2(k,kk)=(z-xopt)'*S*(z-xopt);
    end
end
figure(1); clf reset;
contour(x,y,J',12);
hold on;
plot(xopt(1),xopt(2),'o','markersize',10)
grid

plot([.5,.5],[-8,2],'m--','linewidth',2);
plot([-3,4],[2,-5],'m--','linewidth',2);
plot([-1,4],[2,-6],'m--','linewidth',2);


%figure(2)
%mesh(x,y,J');