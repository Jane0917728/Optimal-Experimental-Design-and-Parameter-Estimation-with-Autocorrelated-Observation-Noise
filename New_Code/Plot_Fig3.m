% Drawing the posterior for the correlated noise, IID noise and
% misspecified noise (assummed IID but it's correlated noise)
r=0.2; K=50; C0=4.5;

SigmaC=3;% The standard deviation of IID noise
 
t_initial=0;
TF=80+t_initial;
set_numpoints=[3,4,5,6,7,8,9,10];
l_set=length(set_numpoints);
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  
symbols = ['o', 'p', '+', '^', 'x', 's', 'd',  'v', '<', '>', 'h'];  

%Estimate the parameters 'prm_interest' from 'prm_name_all'.
prm_all=[r K C0]'; prm_name_all = {'r', 'K', 'C_0'};
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);
lb=[0.1,30,0.01];
ub=[0.4,60,15]; 
numpts=70; 
t_index=3; % indicate choosing 5 measurement points
 
 
est_matrix=zeros(3,3);
 

 
type={'Evenly','FIM','Sobol'};
types=length(type);
 
Z=randn(1,set_numpoints(t_index));

 IID_noise=Z(1:end)*SigmaC;
 initial=[0.4,70,20]; 
 SIG_IID=eye(set_numpoints(t_index),set_numpoints(t_index))*SigmaC^2;  

 for k=1:types    
      fixed1=[0,0,0]; fixed=fixed1;
      SIG=SIG_IID;
      noise=IID_noise;
      if k==1
        %evenly          
         for i=1:l_set
             n_s=set_numpoints(i);
             t_best{i}=linspace(t_initial,TF,n_s);
         end
         
      end
      
      if k==2         
          %IID FIM
            t_best={[6.34993960242455,16.7861647871251,79.9989707801680]
            [6.20976361587004,15.2750515082241,18.0693101564831,79.9952271917983]
            [6.27576505071692,15.7021142959580,17.7021201416564,77.9958861320943,79.9984837351004]
            [5.34128918048555,7.34129105358165,15.8462836262108,17.8462859513070,77.9976752248485,79.9992990846552]
            [5.33730202955380,7.33730312219993,15.8329472137136,17.8329484647109,75.9947523674362,77.9972977774983,79.9988672374438]
            [5.19738604749488,7.19738822874342,14.7063467920120,16.7063514595878,18.7063549261692,75.9830308811366,77.9968807504694,79.9995856922277]
            [4.34607385051873,6.34608727246288,8.34610274025137,14.9428040320692,16.9428216757022,18.9428391065568,75.9159964586665,77.9567198725748,79.9810384182459]
            [4.33685018868102,6.33685084888361,8.33685160332384,14.9132176418497,16.9132184703779,18.9132192649014,73.9951715833340,75.9972308260152,77.9983739074229,79.9992444413468]
            }; 
        
     
      end
      if k==3           
           %IID Sob               
                    t_best={[8,16,80]
            [8,10,16,80]
            [8,10,16,74,80]
            [8,10,16,18,74,80]
            [8,10,16,18,72,74,80]
            [8,10,14,16,18,72,74,80]
            [6,8,10,14,16,18,72,74,80]
            [6,8,10,14,16,18,72,74,78,80]
            };
            
      end
        t=t_best{t_index};
        N_real=Nfunction(r,K,t,C0);
        measurements=N_real+noise;
        [mlevalue,max_likelihood] =fmincon_likenew(fixed,initial,lb,ub,SIG,measurements,t);

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
                [minimizer,max_likelihood,~,~] =  fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t);

                %     [minimizer,MSE_v(param,i),max_ls(param,i),~,~,~,~] =fmincon_likelihood2(initial_likelihood,fixed_params,lb,ub,N0,t,NData,Sigma);
                minimizers(fixed==0,i)=minimizer;
                max_ls(param,i)=-max_likelihood;
        %         mse=(Nfunction(minimizer(1),minimizer(2),t,N0)-Nfunction(r,K,t,N0)).^2/num_steps;
            end

            for i=mle_idx-1:-1:1
        %         fprintf('Optimizing for %s=%.3f\n',param_names{param},param_vals(param,i));
                initial_likelihood(fixed==0)=minimizers(fixed==0,i+1);
                initial_likelihood(param)=param_vals(param,i);
                [minimizer,max_likelihood,~,~] = fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t);
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
         subplot('Position',[0.05+0.043+param*0.225,0.14,0.21,0.75]);%   [left bottom width height]
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
        
      
 
   
