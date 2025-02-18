function [results,t_opt]=Opt_con_fmin_OU(r,K,N0,phi,min_int,max_int,x0,lb,ub)

 % decision varaible's upper boundary 

% x0=[0.100000000000000,14.7985046931648,80.0999999999998];
f1=@(x) obj_conOU(r,K,phi,N0,x);
constrain=@(x) myConstraints(min_int,max_int,x); 
    
options=optimoptions('fmincon','Algorithm','interior-point');
%     options = optimoptions('fmincon', 'Display', 'off');
%     options.Display='iter';
options.Display='off';
options.Diagnostics='off';
options.MaxFunctionEvaluations=6000;
problem.objective=f1;
problem.x0=x0;
problem.solver='fmincon';
problem.lb=lb;
problem.ub=ub;
problem.options=options;
problem.nonlcon = constrain;

[solution,detFIM,exitflag,fmincon_output,~,grad,hessian] = fmincon(problem); 
%  disp(fmincon_output);




results=detFIM;
t_opt=solution; 
end
 

 