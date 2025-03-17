 
 set_numpoints=[3,4,5,6,7,8,9,10];
% Define colors and markers
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  
% symbols = ['o', '+', '*', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'];  
symbols = ['o', 'p', '+', '^', 'x', 's', 'd',  'v', '<', '>', 'h'];  
% x_vec={x_OU_ef,x_OU_FIM,x_IID_ef,x_IID_FIM};
% x_vec={x_OU_Sob,x_OU_FIM,x_IID_Sob,x_IID_FIM};
%  t_vec={t_best,t_best};
 

t_vec={t_best};
% % t_vec={t_best};
t_vec={t_bestSob};
% t_vec={t_best_OUSob};
% t_best=t_best_OUSob;
% t_best=t_best_IIDFIM;
% t_best=t_best_OUFIM;
% t_best={[0.500000009122466,14.3912834658051,31.3565763516039]
% [0.500000009166573,14.4902279979727,29.3619388900281,80.4999958900493]
% [0.500000009479948,11.9015998570552,16.4967492841004,29.7741552281267,80.4999957768247]
% [0.500000009635325,5.91787758314801,12.6660187377760,16.9399994011307,29.8423846550065,80.4999956845064]
% [0.500000038646330,5.82189047644720,12.4937605881981,16.5466886800700,25.5233736374683,34.9336929369821,80.4999841245735]
% [0.500000038811548,5.43777843035908,11.2840634295954,14.5084823910181,18.0045206638338,26.2521339278123,35.6932614189942,80.4999841981693]
% [0.500000009717110,5.36763645483897,11.0469987918951,14.1865924827194,17.4044049679731,23.7738449533160,29.9805182731035,40.6100546533478,80.4999962260149]
% [0.500000007794834,4.91706970987467,9.54765922318163,12.6978431911896,15.2960153254392,18.3734020205038,24.6936160254598,30.8729760280738,42.2372504039319,80.4999970325820]
% };
%   t_vec={t_best};
t1=linspace(t_initial,TF,100)';
Nvalue=Nfunction(r,K,t1,N0); 
% t_best=t_best_OUSob;
N_normalized=Nvalue/Nvalue(end)*(set_numpoints(end)-set_numpoints(1))+set_numpoints(1);
noise_name={'Ornstein-Uhlenbeck','Ornstein-Uhlenbeck','IID','IID'};
% index_name={'eFAST','FIM','eFAST','FIM'};
index_name={'Sobol','FIM','Sobol','FIM'};

fig=figure('Position',[50 50 1350 370],'color','w');
  %   [left bottom width height]
for j=1:length(t_vec)
     subplot('Position',[0.06+(j-1)*0.31,0.13,0.28,0.8]);
     t_best=t_vec{j};
%     x_best=x_vec{j};
%     [num_row,num_time]=size(x_OU_ef);
%     vectors = cell(size(x_best, 1),1); 
%     for i = 1:size(x_best, 1)
%         row = x_best(i, :);  
%         vectors{i} = t(row == 1);  %Store vectors in a cell array for easier iteration
%     end 
  
    hold on; %   Keep the plot from clearing between vector plots
    y_positions =set_numpoints; 
    for i = 1:length(t_best)
        vector = t_best{i};
        y = repmat(y_positions(i), 1, length(vector)); 
        color = colors(mod(i-1, length(colors)) + 1);
        symbol = symbols(mod(i-1, length(symbols)) + 1);   
        plot(vector, y, symbol, 'MarkerSize',8,'LineWidth', 1.2);       
    end
%     plot(t1,N_normalized,'-g','LineWidth',0.8);
    xlim([-1, TF+1]); %  
    ylim([set_numpoints(1)-0.5, set_numpoints(end)+0.5]); %     
    xlabel('Optimized Points on the Time Axis');   
  if j==1
       ylabel('Number of points');
      title('Optimized points using FIM under IID noise');
  end
  if j==2
       title('Optimal points using FIM under OU noise with \phi=0.02');
  end
  if j==3
       title('Optimized points  with Sobol for OU noise with phi=0.02');
  end
%      title(sprintf('Optimized points with %s for %s noise', index_name{j},  noise_name{j}));
%     title('Optimized points with FIM for OU noise');
%  title('Optimal points using FIM under OU noise with \phi=0.02');
%   title('Optimized points  with Sobol for OU noise with phi=0.02');
%   title('Optimized points using FIM under IID noise with \Delta_t=0.5');
%   title('Optimized points using FIM under IID noise with \Delta_t=2');
end

 
 


