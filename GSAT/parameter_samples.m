function [ps,com_ps]= parameter_samples(num_samples,lb,ub,num_params)

s = sobolset(num_params, 'Skip', 1e3, 'Leap', 1e2);
s_com = sobolset(num_params, 'Skip', 1e2, 'Leap', 1e3);
sobol_ps = net(s, num_samples); 
sobol_com = net(s_com, num_samples); 

ps=zeros(num_samples,num_params);
com_ps=zeros(num_samples,num_params);

    for param=1:num_params
    
        ps(:,param)=lb(param)+(ub(param)-lb(param)).*sobol_ps(:,param);
        com_ps(:,param)=lb(param)+(ub(param)-lb(param)).*sobol_com(:,param);
    
    end

end