function [minimizer,max_likelihood,exitflag,fmincon_output] = fmincon_likenew(fixed,initial,lb,ub,SIG,measurements,t_opt)
param_str='[';
paramcount=1;

for i=1:size(fixed,2)
 if ~fixed(i)
            param_str=strcat(param_str,'x(',num2str(paramcount),'),');
            paramcount = paramcount+1;
  else
            param_str=strcat(param_str,num2str(initial(i)),',');
  end
end

param_str=strcat(param_str,']'); 
% L=est_likelihood_fix(param_list,N0,measurements,num_time,t_index,t1,Sigma,rho,phi)
% f_str=strcat('f=@(x) est_likelihood_fix(',param_str,',N0,mu,SIG,measurements,t1,t_index);');
f_str=strcat('f=@(x) est_likenew(',param_str,',SIG,measurements,t_opt);');
% est_likelihood(param_list,N0,mu,SIG,measurements,t_opt)
% f_str=strcat('f=@(x) est_likelihood_fix(',param_str,',measurements,num_time,t_index,t1,Sigma,rho,phi);');
eval(f_str);
options=optimoptions('fmincon','Algorithm','interior-point');
%     options = optimoptions('fmincon', 'Display', 'off');
%     options.Display='iter';
    options.Display='off';
    options.Diagnostics='off';
    options.MaxFunctionEvaluations=2000;
    problem.objective=f;
    problem.x0=initial(fixed==0);
    problem.solver='fmincon';
    problem.lb=lb(fixed==0);
    problem.ub=ub(fixed==0);
    problem.options=options;

    [minimizer,max_likelihood,exitflag,fmincon_output,~,grad,hessian] = fmincon(problem);


  
end
