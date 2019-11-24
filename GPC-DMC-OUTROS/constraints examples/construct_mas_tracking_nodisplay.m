%  SWITCH OFF SCREEN OUTPUT TO SPEED UP ALGORITHM
%
% CONSTRUCT_MAS  (Modified from original by Bert Pluymers)
%
% constructs the maximal output admissable set S for a linear system 
% x(k+1)= Phi x(k) and subject to sample constraints 
% defined by A_y1 * x <= b_y1   and A_y2 * x <= b_y2
% The former of these does not change under the dynamics Phi.
%
%    Garbage_collection :
%
%       0 = do not remove redundant constraints
%       1 = remove redundant constraints after constructing the MAS
%       2 = remove redundant constraints during and after constructint the MAS
%          (speeds up the algorithm in some cases)
% 
% Return values :
%
%    A_S, b_s : matrices defining the constructed invariant set S 
%                   S = {x | A_S * x <= b_S} 
%
function varargout = construct_mas_tracking(Phi,A_y1,b_y1,A_y2,b_y2,garbage_collection)
  
  small = 10^-11;
  
  options = optimset('Display','off');
  
  % generate initial guess for MAS (sample constraints)
  A_S     = [A_y1;A_y2];
  b_S     = [b_y1;b_y2];
  depth_S = zeros(size(b_S));
%  fprintf('Intialization, #constraints =%4.0f\n',size(A_S,1));
  
  % go over the rows of A_S matrix, and add constraints when necessary until
  % it is invariant 
  row=length(b_y1)+1;  %%% begin with 1st row of A_y2
  nf1=row;
  iter=0;
  while row <= size(A_S,1)
    iter=iter+1;

      [x, fval,exitflag] = linprog(-A_S(row,:)*Phi,A_S,b_S,[],[],[],[],[],options);
      if exitflag==-3;
          %disp('Possible problem with admissible set - solution is unbounded');
         % disp(['size of A_S is currently ',num2str(length(b_S))])
          fval=-b_S(row,:)-10;
      end
      if fval<-b_S(row,:)-small
        %fprintf('Constraint added, #constraints =%4.0f\n',size(A_S,1)+1);
        A_S     = [A_S; A_S(row,:)*Phi];
        b_S     = [b_S; b_S(row,:)];
        depth_S = [depth_S; depth_S(row,:)+1];
      end
    
    row = row + 1;
    if (garbage_collection>1) & (mod(iter,10)==0)
      for j=size(A_S,1):-1:nf1;  %%% Must not remove G1*x-f1
        [x, fval,exitflag] = linprog(-A_S(j,:),A_S([1:j-1 j+1:end],:),b_S([1:j-1 j+1:end],:),[],[],[],[],[],options);
        if and(fval>-b_S(j,:)-small,exitflag==1)
         % fprintf('Constraint %4.0f removed, #constraints =%4.0f\n',j,size(A_S,1)-1);
          A_S     = A_S([1:j-1 j+1:end],:);
          b_S     = b_S([1:j-1 j+1:end],:);
          depth_S = depth_S([1:j-1 j+1:end],:);
          if j<row
            row=row-1;
          end
        end
      end
    end
  
  end
  
  
  if garbage_collection>0
   % fprintf('Cleaning up\n');
    % clean up invariant set (remove redundant constraints)
    for i=size(A_S,1):-1:1;
      [x, fval,exitflag] = linprog(-A_S(i,:),A_S([1:i-1 i+1:end],:),b_S([1:i-1 i+1:end],:),[],[],[],[],[],options);
      if and(fval>-b_S(i,:)-small,exitflag==1)
     %   fprintf('Constraint %4.0f removed, #constraints =%4.0f\n',i,size(A_S,1)-1);
        A_S     = A_S([1:i-1 i+1:end],:);
        b_S     = b_S([1:i-1 i+1:end],:);
        depth_S = depth_S([1:i-1 i+1:end],:);        
      end
    end
  end
  
  varargout{1} = A_S;
  varargout{2} = b_S;
  if nargout>2
    varargout{3} = depth_S;
  end
