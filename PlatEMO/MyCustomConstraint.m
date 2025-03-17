function [c,ceq]=MyCustomConstraint(x, nOnes)
 c=[sum(x(1:end-2))-nOnes;-sum(x(1:end-2))+nOnes];
% ceq=sum(x)-nOnes;
 %c=sum(x)-nOnes;
ceq=[];
end