%Main program for optimization of the number of the observed data, the
%series t and N0.tf FIM local index
% clc;
% clear; 
%initialization
%not optimize N0 nor tf
r=0.2;
K=50;
num_time=20;
num_points=4;
% tf=60;N0=1;
% fN=@(r,K,t,N0)  K*N0./((K(1+uK(t))-N0).*exp(-r*(1+ur(t)).*t)+N0);
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0);
rho=0;
phi=0;
% rho=0.98;
% phi=0.98;

%%%Generate ARMA(1,1) noise 

% measurements=f(r,t,K,N0)+epsilon; 



 
% x_initial=ones(num_time,1);
% nvars=length(x_initial);
% nvars=num_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower and upper bounds
lb=[zeros(num_time,1)',1,10];
ub=[ones(num_time,1)',20,80];
% intcon=1:num_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  f1=@(x) obj1(r,K,num_time,rho,phi,x); 
 
 f1=@(x) obj1_tfN01(r,K,num_time,rho,phi,x);
 f2=@(x) obj2_tfN01(r,K,num_time,rho,phi,x);

%constraint
g1=@(x) sum(x(1:num_time))-num_points;
g2=@(x) -sum(x(1:num_time))+num_points;
%encoding
enc_b = zeros(1,num_time)+2;
% enc_c = [1 1];
% enc = [enc_b enc_c];
enc=[enc_b,1,1];

%platemo('save',1,'algorithm',@SparseEA,'objFcn',{f},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub)
%   [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
 %[Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%   [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1,f2},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
  [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1,f2},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%   [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% nonlcon = @(x) MyCustomConstraint3(x, num_points); % constraint that enforce the selection of exactly 'n' ones. 
 
 

% Sovle this mixed integer programming by intlinprog
% f=@(x) obj_autocorrelated3(r,K,num_time,rho,phi,N0,tf,x); 
%  options = gaoptimset('EliteCount', 5);

%  options = gaoptimset('MaxGenerations', 100, 'StallGenLimit', 20);

%  options = optimoptions('ga', 'PopulationSize', 100, 'MaxGenerations', 100, 'CrossoverFraction', 0.8, 'MutationRate', 0.01, 'Display', 'iter');
%   options = optimoptions('ga','ConstraintTolerance',1e-6,'PopulationSize',200,'MaxGenerations',5000, 'StallGenLimit', 50, 'CrossoverFraction', 0.7,'EliteCount', 5);
% %   options = optimoptions('ga','ConstraintTolerance',1e-6,'Display','iter');
%   [x, fval, exitFlag, output] = ga(f, nvars, [], [], [], [], lb, ub, nonlcon, intcon,options);
%   [x, fval, exitFlag, output] = ga(f, nvars, [], [], [], [], lb, ub, nonlcon, intcon);
%  disp(['output',output.message]);
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
% disp(fval);
 % full display
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dicision parameters from x

% tf=x(end)
% x=Dec(end,:);
[obj_smalleast, obj_index] = min(0.5*Obj(:,1)/max(Obj(:,1))+0.5*Obj(:,2)/max(Obj(:,2)));
x=Dec(obj_index,:);
obj_com=obj_smalleast
obj11=obj1_tfN01(r,K,num_time,rho,phi,x)
obj21=obj2_tfN01(r,K,num_time,rho,phi,x)
variables=x(1:num_time)';
t_index=find(variables==1);

tf=x(end)
N0=x(end-1)
 t1=linspace(0,tf,num_time);
t=t1(t_index); 

 
N_selection=fN(r,K,t,N0);
N=fN(r,K,t1,N0);
figure;
plot(t1, N, '-',t,N_selection,'o','LineWidth', 2);  
xlabel('time');
ylabel('Density');
legend('Model curve', 'Selected points');
% title('IID noise');
title('ARMA(1,1) noise with \rho=0.98, \phi=0.98')

% points_num=sum(variables)
% result=[t;fval];
% save("D:\Matlab\Autocorrelated_optimize\data\optimized_ARMA.mat","t",'-append')
% filename='D:\Matlab\Autocorrelated_optimize\data\time_in.mat';
% filename = fullfile('D:\Matlab\Autocorrelated_optimize\data\', 'matlab.mat');
% save('matlab.mat', 'result', '-append');

 