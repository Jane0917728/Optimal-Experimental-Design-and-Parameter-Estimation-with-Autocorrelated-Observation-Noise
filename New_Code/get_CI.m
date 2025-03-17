function [CI_r,CI_K,CI_N0,MSE1] = ...
         get_CI(r,K,N0,numpts,SIG,lb,ub,measurements,t_opt,num_params)

%     rho=0; phi=0;
%    num_steps = 3;   % number of observed points
%
%  num_time=length(t1);
%  num_steps=length(t);
%  N=Nfunction(r,K,t,N0);
CI_r=zeros(3,1);
CI_K=zeros(3,1); 
CI_N0=zeros(3,1);
initial=[0.15,60,3];
% initial2=[0.05,800];
% num_params=3;
fixed1=[0,0,0]; fixed=fixed1;
% fmincon_Like_fix(fixed,initial,lb,ub,rho,phi,Sigma,N0,num_time,measurements,t1,t_index)
[mle,max_likelihood] =fmincon_likenew(fixed,initial,lb,ub,SIG,measurements,t_opt);
CI_r(2)=mle(1);
CI_K(2)=mle(2);
CI_N0(2)=mle(3);
Max_l=max_likelihood;
N_est = Nfunction(mle(1,1),mle(1,2),t_opt,mle(1,3));
N_real=Nfunction(r,K,t_opt,N0);
MSE1=sum((N_est-N_real).^2)/length(N_real);


%% profile likelihood
% An array of parameters of interest (assign the value)
param_vals=zeros(num_params,numpts);
% Likelihood array
max_ls=zeros(num_params,numpts); 
% Optimized parameter array
minimizers=zeros(num_params,numpts); 
% MSE_v=zeros(num_params,numpts);
% Add global optimal values
% optimal_param_vals=initial2';
optimal_param_vals=mle';
param_vals=[param_vals,optimal_param_vals];
% param_name=['r','K'];

%%%for r
 initial_likelihood=optimal_param_vals;
for param=1:num_params    
% The parameters are equally divided and the index of the globally optimal parameter is obtained
    param_vals(param,1:numpts)=linspace(lb(param),ub(param),numpts); % assign values to a parameter.
    param_vals(param,:)=sort(param_vals(param,:));
    [~,mle_idx]=min(abs(param_vals(param,:)-optimal_param_vals(param)));% find the index of the optimal parameter
   fixed=fixed1;
   fixed(param)=1;    
   minimizers(fixed==0,mle_idx)=optimal_param_vals(fixed==0);
%    minimizers(:,mle_idx)=[optimal_param_vals(2),optimal_param_vals(1)];
   max_ls(param,mle_idx)=-Max_l;
%    MSE_v(param,mle_idx)=mse;

    for i=mle_idx+1:numpts+1
          % if optimize K, given value of r, let initial K=optimal value of K in previous step 
       initial_likelihood(fixed==0)=minimizers(fixed==0,i-1); 
              % Let r=given value defined in param_value.
        initial_likelihood(param)=param_vals(param,i);
        [minimizer,max_likelihood,~,~] =  fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t_opt);
  
        %     [minimizer,MSE_v(param,i),max_ls(param,i),~,~,~,~] =fmincon_likelihood2(initial_likelihood,fixed_params,lb,ub,N0,t,NData,Sigma);
        minimizers(fixed==0,i)=minimizer;
        max_ls(param,i)=-max_likelihood;
%         mse=(Nfunction(minimizer(1),minimizer(2),t,N0)-Nfunction(r,K,t,N0)).^2/num_steps;
    end

    for i=mle_idx-1:-1:1
%         fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
        initial_likelihood(fixed==0)=minimizers(fixed==0,i+1);
        initial_likelihood(param)=param_vals(param,i);
        [minimizer,max_likelihood,~,~] = fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t_opt);
         minimizers(fixed==0,i)=minimizer;
         max_ls(param,i)=-max_likelihood;
    end

end


%%

%   num_free_params=sum(1-fixed);
% fig=figure('Position',[100 100 1400 400],'color','w');
% free_param_count=1;
% Stores the zero of the Profile Likelihood curve
zs = cell(num_params,1);
% Store The width of the store confidence interval
% conf_interval=nan(num_params,1);

for param=1:num_params    

    xx=param_vals(param,:);
    yy=max_ls(param,:)-max(max_ls(param,:));    
    zs{param}=interp_zero(xx,yy+1.92);
    if size(zs{param},2) == 2        
             if param==1
                    CI_r(1)=zs{param}(1);
                    CI_r(3)=zs{param}(2);
             elseif param==2
                    CI_K(1)=zs{param}(1);
                    CI_K(3)=zs{param}(2);
             else
                    CI_N0(1)=zs{param}(1);
                    CI_N0(3)=zs{param}(2);
             end             
    else

        if size(zs{param},2) == 1
            if(zs{param} > mle(1,param))
                if param==1
                    CI_r(1)=-0.2;
                    CI_r(3)=zs{param}(1);
                elseif param==2
                    CI_K(1)=10;
                    CI_K(3)=zs{param}(1);
                else
                    CI_N0(1)=-1;
                    CI_N0(3)=zs{param}(1);
                end    
            else
                if param==1
                    CI_r(1)=zs{param}(1);
                    CI_r(3)=1;
                elseif param==2
                    CI_K(1)=zs{param}(1);
                    CI_K(3)=110;
                else
                    CI_N0(1)=zs{param}(1);
                    CI_N0(3)=30;
                end    
            end

        else 
            if param==1
                    CI_r(1)=-0.2;
                    CI_r(3)=1;
            elseif param==2
                    CI_K(1)=10;
                    CI_K(3)=110;
            else
                    CI_N0(1)=-1;
                    CI_N0(3)=30;
             end            
        end
    end

end
