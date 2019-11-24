%%% Find a MAS given the following dynamics
%%%     F x <=t
%%% Process is x(k+1)=Ax(k)
%%% Constraints at each sample are C x = [C1;C2]x <=[f1;f2]=f
%%%
%%%  Note constraints in rows C2 do not change so are excluded from
%%%  iterative step.
%%%
%%% Uses a simple while loop and simple method based 
%%% on gradually increasing the number of rows of a matrix inequality
%%% does not remove redundant constraints
%%%  F=[C;C1 A;C1 A^2;...] <=[f;f1;f1;...]
%%%
%%%  [F,t]=findmas_tracking(A,C1,C2,f1,f2)

function [F,t]=findmas_tracking(A,C1,C2,f1,f2)

%%%% initial set for k=0;
F=[C1;C2];
t=[f1;f2];
val=zeros(size(f1));
An=C1;
nc=size(C1,1);  %%% number of constraints in each sample
Inc=eye(nc);

%%% Switch of display
opt = optimset('linprog');
opt.Diagnostics='off';    %%%%% Switches of unwanted MATLAB displays
opt.LargeScale='off';     %%%%% However no warning of infeasibility
opt.Display='off';


cont=1;
while cont==1;
    An=An*A;   %%% forms new block An=C1*A^n
    
    %%% inequalities to maximise - check each row of new constraints
    %%%  e(j)^T An x <= t(j)
    for j=1:nc;
        vec=-Inc(j,:)*An;
        [x,vv,exitflag]=linprog(vec,F,t,[],[],[],[],[],opt);
        if exitflag==1; val(j)=vv;
        elseif exitflag==-3; %%% indicates solution unbounded
            val(j)=-f1(j)*2; %% marks this row as needed
        end
    end
   
    %%% Extra rows are needed if any values of val exceed f
    %%% so add them in
    if any(-val>f1);
        
    F=[F;An];  %% Add extra block to F
    t=[t;f1];   %% Add extra block to t
        else
        cont=0;  %%% finish loop as all new rows redundant
    end
end


    