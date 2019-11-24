%%%% Determine and sketch MAS
%%% Process is x(k+1)=Ax(k)
%%% Constraints are Cx <=f

%%%% SISO - result obvious
A=0.8;
C=[1;-1];
f=[1;0.5]
[F,t]=findmas(A,C,f)

AA=[1;A;A^2;A^3;A^4;A^5;A^6;A^7;A^8];
figure(1); clf reset
plot([0,8],[1 1],'k--',[0,8],[-0.5,-0.5],'k--');
hold on
plot([0:8]',AA*1,'r',[0:8]',AA*(-0.5),'b');
ylim([-0.6,1.2]);
legend('Upper constraint','lower constraint','simulation 1','simulation 2');


