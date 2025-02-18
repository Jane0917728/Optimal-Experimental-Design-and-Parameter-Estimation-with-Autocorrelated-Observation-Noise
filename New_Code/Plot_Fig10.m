% Define the file names of the workspace files
% workspace_files = {'CI_OU_1.1.mat', 'CI_OU_phi0.1sig1.1.mat', 'CI_OU_phi0.15sig1.1.mat', 'CI_OU_phi0.3sig1.1.mat'};
workspace_files = {'CI_opt_OUSig2_phi0.02_500_2.mat', 'CI_opt_OUSig2_phi0.13_500_2.mat', 'CI_opt_OUSig2_phi0.3_500.mat'};

% Initialize an empty array to store the combined CI data
combined_CI = [];

% Loop through each workspace file and load the CI variable
for i = 1:length(workspace_files)
    % Load the workspace file
    data = load(workspace_files{i}, 'CI');
    
    % Check if the CI variable exists in the workspace
    if isfield(data, 'CI')
        % Concatenate the CI data to the combined array
        combined_CI = [combined_CI, data.CI(:,[1:3,7:9])];
    else
        warning('Variable CI not found in %s', workspace_files{i});
    end
end

% Optionally save the combined CI data to a new workspace file
save('combined_CI.mat', 'combined_CI');


type_name = {'IID FIM', 'OU FIM'};
num_types=length(type_name);
phi_value={'0.02','0.13','0.3'};
num_phi=length(workspace_files);

CI_all=combined_CI;
% colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']; 
% low_index = linspace(0, 0, l_set);
% up_index = linspace(0, 0, l_set);
% mid_index = linspace(0, 0, l_set);

   index = 1:l_set;
    low_index  = 4 * index - 3;
    mid_index  = 4 * index - 2;
    up_index  = 4 * index - 1;
 
% para_dis=[0.02,4,1];
% Define color list, can be modified as per your need
colors = lines(num_types); % Using lines colormap to generate distinct colors
for j = 1:num_param
    fig = figure('Position', [20 20 1000 320], 'color', 'w');    
    ii=[1,4,7,10,13,16]; 
    para_max=max(CI_all(3,ii+(j-1)));
     para_min=min(CI_all(1,ii+(j-1)));
     para_err1=[0.01,1.5,0.2];
      para_err2=[0.01,2,0.2];
    for k=1:num_phi
          subplot('Position', [0.04 + (k - 1) * 0.32, 0.14, 0.3, 0.78]);
             hold on;
             ylim([para_min-para_err1(j),para_max+para_err2(j)]);    
             xlim([set_numpoints(1),set_numpoints(end)+0.5]);        
         
        % Plot additional lines on the chart
%         plot(set_numpoints', CI_param(:, (m - 1) * num_param + j), lineStyles{m}, 'MarkerSize', 8, 'LineWidth', 1.5);
          
          title(sprintf('phi=%s', phi_value{k}));  
          middle_valuesIID = CI_all(mid_index,  2*(k - 1) * num_param + j); % Middle values are in the 2nd column
          lower_errorsIID = middle_valuesIID - CI_all(low_index,2* (k - 1) * num_param + j); % Difference between middle and lower bound (column 3)
          upper_errorsIID = CI_all(up_index,2* (k - 1) * num_param + j) - middle_valuesIID; % Difference between upper bound (column 1) and middle

          
          middle_valuesOU= CI_all(mid_index, 2* (k - 1) * num_param + 3+j); % Middle values are in the 2nd column
          lower_errorsOU = middle_valuesOU - CI_all(low_index, 2*(k - 1) * num_param + 3+ j); % Difference between middle and lower bound (column 3)
          upper_errorsOU = CI_all(up_index, 2*(k - 1) * num_param + 3+ j) - middle_valuesOU; % Difference between upper bound (column 1) and middle

        % Plot error bars
%         errorbar(set_numpoints, middle_valuesIID, lower_errorsIID, upper_errorsIID, '*--',...
%            'MarkerFaceColor', 'c',...
%             'LineWidth', 1.5, 'MarkerSize', 8);
%         errorbar(set_numpoints+0.2, middle_valuesOU, lower_errorsOU, upper_errorsOU, 'o--',...
%             'MarkerFaceColor', 'c','LineWidth', 1.5, 'MarkerSize', 8);
          plot([set_numpoints,set_numpoints(end)+0.2], linspace(prm(j), prm(j), length(set_numpoints)+1),...
               'k--', 'MarkerSize', 8, 'LineWidth', 2);

       
         % Plot the second error bar with enhanced style
        
            
            % Plot the first error bar with enhanced style    
            
          errorbar(set_numpoints, middle_valuesIID, lower_errorsIID, upper_errorsIID, '*-',...
                'Color', [0.2, 0.6, 0.8],... % Light blue color for better aesthetics
                'MarkerFaceColor', [0.2, 0.6, 0.8],... % Fill marker with light blue
                'LineWidth', 1.5, 'MarkerSize', 8);

               errorbar(set_numpoints + 0.2, middle_valuesOU, lower_errorsOU, upper_errorsOU, 'o-',...
                'Color', [0.9, 0.4, 0.1],... % Orange color for better contrast
                 'LineWidth', 1.5, 'MarkerSize', 8);
%                 'MarkerFaceColor', [0.9, 0.4, 0.1],... % Fill marker with orange
               

           if k == 1
               ylabel(sprintf('Parameter %s', prm_name{j}));
           end
           if k == 1                
             legend('real','IID FIM', 'OU FIM', 'Location', 'NorthEast'); % Add legend
 
           end
        xlabel('Number of points');
        % Ensure all x-axis points are displayed
        xticks(set_numpoints);
          axis('square');  
        grid on;
    end
    hold off;
end
