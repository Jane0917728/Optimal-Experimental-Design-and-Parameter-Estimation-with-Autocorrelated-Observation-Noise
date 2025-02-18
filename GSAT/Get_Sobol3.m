% Function for computing the Sobol's indice
function [T,S,St]=Get_Sobol3(param_list, lb,ub,t,num_MC)
% 
% Citation: Cannavo' F., Sensitivity analysis for volcanic source modeling quality assessment and model selection, Computers & Geosciences, Vol. 44, July 2012, Pages 52-59, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2012.03.008.
 

% N0=1;
% tf=80;
% num_time=20;
% n_step=4;
% % t=linspace(0,tf,num_time)';
% t=20;
% create a new project 
pro = pro_Create();

% add 2 input variables with a pdf uniformely distributed in [0 1]
% pro = pro_AddInput(pro, @(N)pdf_Uniform(N, [0 1]), 'param1');
% pro = pro_AddInput(pro, @(N)pdf_Uniform(N, [0 1]), 'param2');
 

% add to the project 2 input variables, named param*, distributed in the 
% range [0 1] and indicate that the variables will be sampled following a 
% Sobol set 
num_param=length(param_list);
for i=1:num_param
    pro = pro_AddInput(pro, @()pdf_Sobol([lb(i) ub(i)]), param_list(i));
%     pro = pro_AddInput(pro, @()pdf_Sobol([lb(2) ub(2)]), 'K');
end



% set the model, and name it as 'model', to the project 
% the model is well-known as "Sobol's function" 
pro = pro_SetModel(pro, @(x)Nfunction3(x,t), 'model');

% calculate the real analytical values of sensitivity coefficients for the 
% model "Sobol's function" (3.0.1).
% [D, Si] = SATestModel([r,K,t,N0]);

% set the number of samples for the quasi-random Monte Carlo simulation
pro.N = num_MC;

% initialize the project by calculating the model at the sample points
pro = GSA_Init(pro);

% calculate the total global sensitivity coefficient for the set of 2 variables
% ('param1 and param5) (see sections 2.3 and 2.4)
% [Stot eStot pro] = GSA_GetTotalSy(pro, {'param1', 'param5'});

% calculate the global sensitivity coefficient for the set of all the input
% variables and verify that equals 1
% [S, eS, pro] = GSA_GetSy(pro, {1,2});

% calculate the first order global sensitivity coefficient for the second variable and
% verify that S2 equals the real value Si(2)
for i=1:num_param
     [S(i,1),~,~]=GSA_GetSy(pro, {i});
%      [St(i,1),~,~]=GSA_GetSy(pro, {i});
     [St(i,1),~,~]=GSA_GetTotalSy(pro,{i});
end
T=pro.GSA.D;
% S(i,1)=pro.GSA.D*S(i,1);
% St(i,1)=pro.GSA.D*St(i,1);
% [Sr, eS1, pro] = GSA_GetSy(pro, {'r'});
% [SK, eS2, pro] = GSA_GetSy(pro, {'K'});
%  [Srtot eStot pro] = GSA_GetTotalSy(pro,{'r'});
%  [SKtot eStot pro] = GSA_GetTotalSy(pro,{'K'});

% calculate the first order global sensitivity coefficients by using FAST
% algorithm and verify that all coefficients equal the real ones in Si
% Sfast = GSA_FAST_GetSi(pro);



end

