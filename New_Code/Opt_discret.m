%Optimizing the time location of the observed data
 function [results,x]=Opt_discret(n_s,Sti,num_time,SIG)
%   Sti= squeeze(Sti(:,:,7));
%   n_s=9;
% %    Sti=S_FIM;
% num_time=41;
%   SIG=SIG_OU;
 f1=@(x)  obj_discret(SIG,Sti,x);
%  f1=@(x)  obj1_Brsimple(dr,dK,SIG_OU,FIM_full_BFIM,x); 
%  f2 =@(x) obj2_Brsimple(dr,dK,SIG_OU,x);

%constraint
g1=@(x) sum(x(1:num_time))-n_s;
g2=@(x) -sum(x(1:num_time))+n_s;
%encoding
enc_b = zeros(1,num_time)+2;% encoding "2" denotes the variables are integers.
% enc_c = [1 1];
% enc = [enc_b enc_c];
enc=enc_b;

lb=zeros(num_time,1)'; %decision varaible's lower boundary
ub=ones(num_time,1)'; %decision varaible's upper boundary, which defines 0 or 1 integer varaibles.

%platemo('save',1,'algorithm',@SparseEA,'objFcn',{f},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub)
%      [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%       [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%       [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1,f2},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%      [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%    [Dec,Obj,Con] = platemo('algorithm',@NSGAII,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%  
% CMA-ES
%     [Dec,Obj,Con] = platemo('algorithm',@NSGAII,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
    [Dec,Obj,Con] = platemo('algorithm',@CSO,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);

% Find the smallest objective value and its index
   [obj_smalleast, obj_index] = min(Obj);
   minIndices = find(Obj == obj_smalleast);

    % If one cannot get an idx satisfying the condition, index is null
    idx = -1;
    %If the element of Con==0, return the idx
    for i = 1:length(minIndices)
        if all(Con(minIndices(i), :) == 0)
            idx = minIndices(i);
            break;
        end
    end
   
    if (idx==-1)
        disp('There is no solution satisfies the constraints');
        x=Dec(1,:);
        results=Inf;
    else
        x=Dec(idx,:);
        results=Obj(idx);
 % [obj_smalleast, obj_index] =min(0.5*Obj(:,1)/max(abs(Obj(:,1)))+0.5*Obj(:,2)/max(abs(Obj(:,2))));
    end

 
   end

 
 