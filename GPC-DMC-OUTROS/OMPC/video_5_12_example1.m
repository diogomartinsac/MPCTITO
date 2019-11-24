%%%%Feasible target computation
A=[0.8,0.1;-0.2,0.9];
B =[0.1;0.8];
C =[1.9 -1];
D=0;
%%%%% constraints
umin=-1;     %%% umin<u<umax
umax=2;
Kxmax=[1 0.2;-0.1 0.4;-1,-0.2;0.1,-0.4];
xmax=[4;4;0.8;2.5];

%%% determine expected steady-states
Mat=inv([C,zeros(1,1);[A-eye(2),B]]);
Kxr=Mat(1:2,1)
Kur=Mat(3,1)

%%%% Test for limits on target
H=[Kur,0;-Kur,0;0,Kur;0,-Kur;Kxmax*Kxr,zeros(4,1);zeros(4,1),Kxmax*Kxr];
epsilon=0.1;
hh=[umax-epsilon;-umin+epsilon;umax-epsilon;-umin+epsilon;xmax-epsilon;xmax-epsilon];
figure(1); clf reset
P=Polyhedron(H,hh);  %%% from mpt3 toolbox
plot(P,'color','y'); 
xlabel('Upper limit on r-d');
ylabel('Lower limit on r-d')
