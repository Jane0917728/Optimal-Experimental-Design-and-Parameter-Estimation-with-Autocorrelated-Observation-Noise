clear all
close all;

set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

r=0.2; 
K=50; 
C0=4.5;
SigmaC=3;% The standard deviation of IID noise
t_initial=0;
TF=80+t_initial;
set_numpoints=[3,4,5,6,7,8,9,10];
l_set=length(set_numpoints);

%Estimate the parameters 'prm_interest' from 'prm_name_all'.
prm_all=[r K C0]'; 
prm_name_all = {'r', 'K', 'C_0'};
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
        % evenly
        for i=1:l_set
            n_s=set_numpoints(i);
            t_best{i}=linspace(t_initial,TF,n_s);
        end

    end
    if k==2
        % IID FIM
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
        % IID Sob
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
        % the parameters are equally divided and the index of the globally optimal parameter is obtained
        param_vals(param,1:numpts)=linspace(lb(param),ub(param),numpts); % assign values to a parameter.
        param_vals(param,:)=sort(param_vals(param,:));
        [~,mle_idx]=min(abs(param_vals(param,:)-optimal_param_vals(param)));% find the index of the optimal parameter
        fixed=fixed1;
        fixed(param)=1;
        minimizers(fixed==0,mle_idx)=optimal_param_vals(fixed==0);
        max_ls(param,mle_idx)=-Max_l;
        
        for i=mle_idx+1:numpts+1
            % if optimize K, given value of r, let initial K=optimal value of K in previous step
            initial_likelihood(fixed==0)=minimizers(fixed==0,i-1);
            % Let r=given value defined in param_value.
            initial_likelihood(param)=param_vals(param,i);
            [minimizer,max_likelihood,~,~] =  fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t);
            minimizers(fixed==0,i)=minimizer;
            max_ls(param,i)=-max_likelihood;
        end

        for i=mle_idx-1:-1:1
            initial_likelihood(fixed==0)=minimizers(fixed==0,i+1);
            initial_likelihood(param)=param_vals(param,i);
            [minimizer,max_likelihood,~,~] = fmincon_likenew(fixed,initial_likelihood,lb,ub,SIG,measurements,t);
            minimizers(fixed==0,i)=minimizer;
            max_ls(param,i)=-max_likelihood;
        end
    end 


    %%
    % Stores the zero of the Profile Likelihood curve
    free_param_count=1;
    % Stores the zero of the Profile Likelihood curve
    zs = cell(num_param,1);
    % Store The width of the store confidence interval
    conf_interval=nan(num_param,1);

    fig=figure('Position',[20 20 1250 320],'color','w');
    subplot('Position',[0.05,0.14,0.23,0.75]);
    hold on;
    y_positions =set_numpoints;
    for i = 1:length(t_best)
        vector = t_best{i};
        y = repmat(y_positions(i), 1, length(vector));
        plot(vector,y,'x','MarkerSize',6,'LineWidth',2);
    end
    xlim([-1,TF+1]); 
    ylim([set_numpoints(1)-0.5,set_numpoints(end)+0.5]);
    xlabel('Observation times');
    ylabel('Number of observation points');
    xticks([0 20 40 60 80]);
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
        subplot('Position',[0.05+0.055+param*0.225,0.14,0.21,0.75]);
        hold on;
        xx=param_vals(param,:);
        yy=max_ls(param,:)-max(max_ls(param,:));
        plot(xx,yy,'-','color',[0,0,0],'LineWidth',2);
        xline(mlevalue(1,param),'--','color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
        plot([min(param_vals(param,:)),max(param_vals(param,:))],[-1.92,-1.92],'color',[0.5, 0.5, 0.5],'LineWidth',2);
        xlabel(prm_name{param});
        ytickformat('%.1f');
        if param==1
            ylabel('Normalized likelihood');
            
        end
        axis('square');
        xlim([min(param_vals(param,:)),max(param_vals(param,:))]);
        ylim([-2.5,0]);
        if param==1
            xtickformat('%.1f');
            xlabel('$r$','Interpreter', 'latex')
        end
        if param==2
            xticks([35 45 55 65]);
            xlabel('$K$','Interpreter', 'latex')
        end
        if param==3
            xlim([0 15])
            xticks([0 5 10 15])
            xlabel('$C_0$','Interpreter', 'latex')
        end
        hold off;

        zs{param}=interp_zero(xx,yy+1.92);
        if size(zs{param},2) == 2
            conf_interval(param)=zs{param}(2)-zs{param}(1);
            fprintf('95%% Confidence interval for param %s is: (intercept at -1.92)\n',prm_name{param});
            fprintf('width=%.4f: [%.4f,%.4f]\n',conf_interval(param),zs{param}(1),zs{param}(2));
        else
            fprintf('Do not have 2 intercepts for param %s, they are:\n',prm_name{param});
            disp(zs{param});
        end      
    end
    print(sprintf('Fig3_%i',k),'-dpng','-r300')
end





