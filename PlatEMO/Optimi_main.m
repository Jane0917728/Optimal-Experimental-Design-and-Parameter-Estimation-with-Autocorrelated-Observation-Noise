%Genetic Algorithm Main program for optimization of the number of the observed data, the series t and N0.
% clc;
clear; 
%initialization
r=0.3;
K=1000;
num_time=200;
num_points=12;
% tf=150;
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0);
rho=0.7;
phi=0.7;
%%%Generate ARMA(1,1) noise 
% Sigma=20;
% v=randn(num_time,1)*20; %column series
% 
% epsilon=zeros(num_time,1);
% epsilon(1)=v(1);
% for tt=2:num_time
%     epsilon(tt,:)=rho*epsilon(tt-1,:)+v(tt,:)+phi*v(tt-1,:);
% end
 
% measurements=f(r,t,K,N0)+epsilon; 


nmax=num_time;
nmin=10;
N0_initial=30;
x_initial=[zeros(num_time,1);N0_initial;100];
nvars=length(x_initial);
%lower and upper bounds
lb=[zeros(num_time,1);2;50]';
ub=[ones(num_time,1);100;150]';
intcon=1:num_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonlcon = @(x) MyCustomConstraint(x, num_points); % constraint that enforce the selection of exactly 'n' ones. 
% Sovle this mixed integer programming by intlinprog
f=@(x) obj_autocorrelated(r,K,num_time,rho,phi,x); 
% options = gaoptimset('MaxGenerations', 100, 'StallGenLimit', 20);
%  options = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 100, 'CrossoverFraction', 0.8, 'MutationRate', 0.01, 'Display', 'iter'，'MutationFcn', {@mutationadaptfeasible, 0.05});
 options = optimoptions('ga','ConstraintTolerance',1e-6,'PopulationSize', 100,'MaxGenerations', 1000, 'StallGenLimit', 20, 'CrossoverFraction', 0.8,'EliteCount', 5);
 [x, fval] = ga(f, nvars, [], [], [], [], lb, ub, nonlcon, intcon,options);
 
%  problem=struct;
%  problem.solver='intlinprog';
%  problem.intcon=intcon;
%  problem.options=options;
%  %problem.f=f(r,K,w1,tf,nmax,x);
%   problem.objective=f;
%  problem.Aineq=A;
%  problem.Bineq=B;
%  problem.intcon=[1:nmax];
%  problem.lb=lb;
%  problem.ub=ub;
%  [x,fval,exitflag,output] = intlinprog(problem);
 
% fprintf('exitflag=%d\n', exitflag); 
 
% disp(x(x~=0));
disp(fval);
 % full display
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dicision parameters from x

%pic
% Delta_t2=x(nmax+1:end-1);
N0=x(end-1);
tf=x(end)
variables=x(1:num_time)';

t1=linspace(0,tf,num_time)';
t2=t1.*variables;
 

% 使用逻辑索引来去掉0元素
t=t2(t2 ~= 0);
N_selection=fN(r,K,t,N0);
N=fN(r,K,t1,N0);
figure;
plot(t1, N, '-',t,N_selection,'o','LineWidth', 2);  
legend('model curve', 'Selection point');
points_num=sum(variables)
num_points
