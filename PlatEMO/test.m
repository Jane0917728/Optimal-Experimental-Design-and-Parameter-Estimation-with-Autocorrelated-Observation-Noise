%Genetic Algorithm Main program for optimization of the number of the observed data, the series t and N0.
% clc;
clear; 
%initialization
r=0.3;
K=1000;
num_time=100;
num_points=12;
nOnes=3;
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

%nonlcon = @(x) MyCustomConstraint(x, num_points); % constraint that enforce the selection of exactly 'n' ones. 
% Sovle this mixed integer programming by intlinprog

%obj
f=@(x) obj_autocorrelated(r,K,num_time,rho,phi,x); 
%constraint
g1=@(x)sum(x(1:end-2))-nOnes;
g2=@(x)-sum(x(1:end-2))+nOnes;
%encoding
enc_b = zeros(1,num_time)+2
enc_c = [1 1]
enc = [enc_b enc_c]


%platemo('save',1,'algorithm',@SparseEA,'objFcn',{f},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub)
[Dec,Obj,Con] = platemo('algorithm',@GA,'objFcn',{f},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub)
%PRO = UserProblem('save',1,'objFcn',f,'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%ALG2 = ABC();
%ALG2.Solve(PRO)
