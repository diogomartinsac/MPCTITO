%SISO example
Dumax=1;
umax=3;
umin=-2;
nu=3;
sizeu=1; %%% SISO case
[CC,dd,du]  = mpc_constraints(Dumax,umax,umin,sizeu,nu)

[CC,dd,du]  = mpc_constraints(Dumax,umax,umin,sizeu,4)

%%% MIMO example
nu=2;
Dumax=[1;0.5];
umax=[1.3;0.9];
umin=[-2;-0.5];
sizeu=2;
[CC,dd,du]  = mpc_constraints(Dumax,umax,umin,sizeu,nu)
