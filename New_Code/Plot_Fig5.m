% Drawing the posterior for the correlated noise, IID noise and
% misspecified noise (assummed IID but it's correlated noise)
r=0.2; K=50; C0=4.5;

SigmaC=3;
phi=0.02;
Sigma=sqrt(SigmaC^2*2*phi);

% Sigma=sqrt(SigmaC^2*(2*phi));
 
t_initial=0;
TF=80+t_initial;
set_numpoints=[3,4,5,6,7,8,9,10];
l_set=length(set_numpoints);
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  
symbols = ['o', 'p', '+', '^', 'x', 's', 'd',  'v', '<', '>', 'h'];  
t_index=3; %Choose 5 measurement points
%Estimate the parameters 'prm_interest' from 'prm_name_all'.
prm_all=[r K C0]'; prm_name_all = {'r', 'K', 'C_0'};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);
lb=[0.05,30,0.01];
ub=[0.4,70,20]; 
numpts=70; 

 %number of measurement points
% t=linspace(t_initial,TF,n_s);
 
est_matrix=zeros(3,3); 
type={'Evenly','FIM','Sobol'};
types=length(type);

 
 initial=[0.3,60,8]; 
 nn_s=set_numpoints(t_index);
    Z=randn(1,nn_s);

 for k=1:types    
      fixed1=[0,0,0]; fixed=fixed1;
     if k==1
        %evenly          
         for i=1:l_set
             n_s=set_numpoints(i);
             t_best{i}=linspace(t_initial,TF,n_s);
         end
         t=t_best{t_index};
         N_real=Nfunction(r,K,t,C0);
        SIG_OU=zeros(nn_s,nn_s); % The covariance matrix 
        for s = 1:nn_s
             for tt = 1:nn_s
        %                 SIG(s, tt) = Sigma^2*(exp(-phi*abs(t_opt(tt)-t_opt(s)))-exp(-phi*(t_opt(tt)+t_opt(s))))/(2*phi);   
                  SIG_OU(s, tt) = Sigma^2*exp(-phi*abs(t(tt)-t(s)))/(2*phi);   
             end   
        end 
       OU_noise=Generate_OU_seeds(phi,Sigma,t,Z);
      end
      
      if k==2         
 
        t_best={[7.63437169214578,15.8005539164688,33.3382268030993]
                    [6.62146611958361e-07,8.36271697631842,16.5826947305807,34.7187852405600]
                    [7.82233269762834e-07,8.24622337773052,16.4095759930023,30.0111301473713,79.9999965970322]
                    [5.98379060179048e-07,7.22223382252997,12.9516071826385,18.3837700093128,31.2295789897367,79.9999967291563]
                    [5.20547784694058e-07,6.16916019822338,10.4399886084285,15.2717087008492,20.2232789321642,32.1979388161369,79.9999968049386]
                    [4.99161991552154e-07,5.84341877015355,9.77031078892515,14.2504120429883,18.4214716599558,25.6611157745421,36.2302888976562,79.9999969965829]
                    [4.79089928281056e-07,5.29297844110465,8.75769687599276,12.4151812253002,15.9968306683055,19.9289913469709,27.4532598526502,38.1092671558424,79.9999970810677]
                    [1.85823376975489e-06,4.83236509493882,7.97717565733385,11.0481370021154,14.3638435126200,17.5188072744450,21.6149109732524,28.9235300217553,39.9099304528990,79.9999890575419]
                    };
        t=t_best{t_index};
         N_real=Nfunction(r,K,t,C0);
          SIG_OU=zeros(nn_s,nn_s); % The covariance matrix 
            for s = 1:nn_s
                 for tt = 1:nn_s
            %                 SIG(s, tt) = Sigma^2*(exp(-phi*abs(t_opt(tt)-t_opt(s)))-exp(-phi*(t_opt(tt)+t_opt(s))))/(2*phi);   
                      SIG_OU(s, tt) = Sigma^2*exp(-phi*abs(t(tt)-t(s)))/(2*phi);   
                 end   
            end 
           OU_noise=Generate_OU_seeds(phi,Sigma,t,Z);
       
      end
      if k==3           
%           t_best={[10,16,32]
%                 [0,10,16,32]
%                 [0,10,16,28,80]
%                 [0,8,12,16,28,80]
%                 [0,8,10,14,18,30,80]
%                 [0,8,10,14,18,26,34,80]
%                 [0,4,8,10,14,18,26,34,80]
%                 [0,4,8,10,14,16,18,26,36,80]
%             };
                    t_best={[10,16,30]
            [2,10,16,30]
            [2,10,16,28,80]
            [0,8,12,16,28,80]
            [0,6,10,14,18,28,80]
            [0,6,10,14,18,26,36,80]
            [0,4,8,10,14,18,26,36,80]
            [0,4,8,10,14,16,20,26,34,80]
            };
            t=t_best{t_index};
            N_real=Nfunction(r,K,t,C0);
            SIG_OU=zeros(nn_s,nn_s); % The covariance matrix 
            for s = 1:nn_s
                 for tt = 1:nn_s
            %                 SIG(s, tt) = Sigma^2*(exp(-phi*abs(t_opt(tt)-t_opt(s)))-exp(-phi*(t_opt(tt)+t_opt(s))))/(2*phi);   
                      SIG_OU(s, tt) = Sigma^2*exp(-phi*abs(t(tt)-t(s)))/(2*phi);   
                 end   
            end 
            OU_noise=Generate_OU_seeds(phi,Sigma,t,Z);
         end
      measurements=N_real+OU_noise;
      [mlevalue,max_likelihood] =fmincon_likenew(fixed,initial,lb,ub,SIG_OU,measurements,t);
    %    prm_all=mlevalue;
       est_matrix(k,:)=mlevalue;
      
      Max_l=max_likelihood;
         %--------------------------------------------------------------------------------------------------------
      %% profile likelihood
      % An array of parameters of interest (assign the value)
      param_vals=zeros(num_param,numpts);
      % Likelihood array
       max_ls=zeros(num_param,numpts); 
       % Optimized parameter array
       minimizers=zeros(num_param,numpts); 
 
       optimal_param_vals=mlevalue';
       param_vals=[param_vals,optimal_param_vals];
 

        %%%for r
         initial_likelihood=optimal_param_vals;
        for param=1:num_param    
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
                [minimizer,max_likelihood,~,~] =  fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG_OU,measurements,t);

                %     [minimizer,MSE_v(param,i),max_ls(param,i),~,~,~,~] =fmincon_likelihood2(initial_likelihood,fixed_params,lb,ub,N0,t,NData,Sigma);
                minimizers(fixed==0,i)=minimizer;
                max_ls(param,i)=-max_likelihood;
        %         mse=(Nfunction(minimizer(1),minimizer(2),t,N0)-Nfunction(r,K,t,N0)).^2/num_steps;
            end

            for i=mle_idx-1:-1:1
        %         fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
                initial_likelihood(fixed==0)=minimizers(fixed==0,i+1);
                initial_likelihood(param)=param_vals(param,i);
                [minimizer,max_likelihood,~,~] = fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG_OU,measurements,t);
                 minimizers(fixed==0,i)=minimizer;
                 max_ls(param,i)=-max_likelihood;
            end
        end %end of num_param


%%


        % Stores the zero of the Profile Likelihood curve
       free_param_count=1;
      % Stores the zero of the Profile Likelihood curve
      zs = cell(num_param,1);
      % Store The width of the store confidence interval
      conf_interval=nan(num_param,1);

        %   [left bottom width height] 
%           fig=figure('Position',[20 20 1000 320],'color','w');
        fig=figure('Position',[20 20 1250 320],'color','w');
        % subplot(1,num_param+1,1);
           
         subplot('Position',[0.05,0.14,0.23,0.75]);
         hold on; %   Keep the plot from clearing between vector plots
         y_positions =set_numpoints; 
         for i = 1:length(t_best)
                vector = t_best{i};
                y = repmat(y_positions(i), 1, length(vector)); 
                color = colors(mod(i-1, length(colors)) + 1);
                symbol = symbols(mod(i-1, length(symbols)) + 1);   
                plot(vector, y, symbol, 'MarkerSize',6,'LineWidth', 1.2);       
         end
%     plot(t1,N_normalized,'-g','LineWidth',0.8);
    xlim([-1, TF+1]); %  
    ylim([set_numpoints(1)-0.5, set_numpoints(end)+0.5]); %     
    xlabel('Optimized Points on the Time Axis');   
    ylabel('Number of points');
      if k==1      
          title('Evenly distributed points');
      end
      if k==2
           title('Optimal points using FIM');
      end
      if k==3
           title('Optimal points using Sobol');
      end
         
         

        for param=1:num_param    

         free_param_count = free_param_count+1;
         subplot('Position',[0.05+0.043+param*0.225,0.14,0.21,0.75]);
  %   [left bottom width height]

            hold on;
            xx=param_vals(param,:);
            yy=max_ls(param,:)-max(max_ls(param,:));
            plot(xx,yy,'-','color',[0,0,0],'LineWidth',2);
            xline(mlevalue(1,param),'color',[1.0 0 0],'LineWidth',1)
            plot([min(param_vals(param,:)),max(param_vals(param,:))],[-1.92,-1.92]);    
            % plot([min(param_vals(param,:)),max(param_vals(param,:))]);
            xlabel(prm_name{param});
            %     ylabel('log(L)');
            if param==1
               ylabel('Normalized Likelihood');
            end
           axis('square');
           xlim([min(param_vals(param,:)),max(param_vals(param,:))]);
           ylim([-2.5,0]);
        %       ylim([-10,0]);
            hold off;    
            zs{param}=interp_zero(xx,yy+1.92);
            if size(zs{param},2) == 2
                conf_interval(param)=zs{param}(2)-zs{param}(1);
                fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',prm_name{param});
                fprintf('width=%.4f: [%.4f,%.4f]\n',conf_interval(param),zs{param}(1),zs{param}(2));
                xline(zs{param}(1), 'color',[0 1 0]);  
                xline(zs{param}(2), 'color',[0 1 1]); 
            else
                fprintf('Do not have 2 intercepts for param %s, they are:\n',prm_name{param});
                disp(zs{param});
                if size(zs{param},2) == 1
                    if(zs{param} > mlevalue(1,param))
                        xline(zs{param}(1), 'color',[0 1 1]);  
                    else
                        xline(zs{param}(1), 'color',[0 1 0]);  
                    end
                end
            end
        end
              %--------------------------------------------------------------------------------------------------------------------------

             
  end
        
      
 
   
