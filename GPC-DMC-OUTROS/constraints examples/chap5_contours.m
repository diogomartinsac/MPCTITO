%%% Contours for a performance index
%%%  J=x'Sx+x'b+c
b=[1 5]';
c=3;
S=[1.2 -0.2;-0.2 0.8];
x=linspace(-6,4,300);
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
title('Contours of fixed J','fontsize',18);

plot([1.5,1.5],[-8,2],'m--','linewidth',2);
plot([-6,4],[7.5,-2.5],'m--','linewidth',2);
plot([-6,4],[14,-15],'m--','linewidth',2);


%figure(2)
%mesh(x,y,J');