  function [c,ceq] = myConstraints(min_int,max_int,t)
     
    s = length(t);     
    % Initialization of inequalities 
    c1 = zeros(s-1, 1);
    c2 = zeros(s-1, 1);
    % t_1 >= 0  transform to  t_1 - 0 <= 0 
    

    % t_i - t_{i+1} + min_int <= 0 and  t_{i+1} - t_i - max_int <= 0
    for i = 1:s-1
        c1(i) =- (t(i+1) - t(i) - min_int);
         c2(i)=t(i+1) - t(i) - max_int;            
    end
%     c3=t(end)-TF;
% c=[c1;c2;c3];
c=[c1;c2];
 ceq=[];

 end
