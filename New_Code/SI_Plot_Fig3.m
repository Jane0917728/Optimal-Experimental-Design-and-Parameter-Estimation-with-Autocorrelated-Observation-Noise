

workspace_files = {'CI_opt_OUSig2_phi0.02_500_2.mat', 'CI_opt_OUSig2_phi0.13_500_2.mat', 'CI_opt_OUSig2_phi0.3_500.mat'};

% Initialize an empty array to store the combined CI data
combined_CI = [];
r=0.2; K=50; N0=4.5;
prm_all=[r K N0]'; prm_name_all = {'r', 'K', 'C_0'};   
prm_interest={'r', 'K', 'C_0'};
prm_idx = find(ismember(prm_name_all, prm_interest));
prm=prm_all(prm_idx);
prm_name=prm_name_all(prm_idx);
num_param=length(prm);

% Loop through each workspace file and load the CI variable
for i = 1:length(workspace_files)
    % Load the workspace file
        data = load(workspace_files{i}, 'CI');
      % Check if the CI variable exists in the workspace
        if isfield(data, 'CI')
            % Concatenate the CI data to the combined array
            CI=data.CI;
            num_param=3;
            set_numpoints=[3,4,5,6,7,8,9,10];
        else
            warning('Variable CI not found in %s', workspace_files{i});
        end

        fig=figure('Position',[20 20 1000 320],'color','w');
          %   [left bottom width height]
          index=1:8;
          index2=[4,8,12,16,20,24,28,32];
          CI_w=CI(index2,:);
        for j=1:num_param
              subplot('Position',[0.06+(j-1)*0.31,0.13,0.265,0.75]);
        %     subplot(1,num_param,j,'Position',[0.05+(j-1)*0.3,0.12,0.26,0.75]);
        %       subplot(1,num_param,j);
        %         plot(set_numpoints, CI_w(kk,:,j+3), 'o--', set_numpoints, CI_w(kk,:,j), 'd--',  'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);
        %     plot(set_numpoints, suc(:,j), '*-',set_numpoints, suc(:,j+num_param*max(0,(num_types-2))), 'o--',...
        %         set_numpoints, suc(:,j+num_param*(num_types-1)), 'd-',  'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);
             plot(set_numpoints(index), CI_w(index,j), '*--',set_numpoints(index), CI_w(index,j+num_param*1), 'x:',...
             set_numpoints(index), CI_w(index,j+num_param*2), 'o--',  set_numpoints(index), CI_w(index,j+num_param*3), 'd:', 'MarkerFaceColor', 'c', 'MarkerSize',8,'LineWidth', 1.5);

              xlim([set_numpoints(index(1)),set_numpoints(index(end))]);
              xlabel('Number of Points');
             if j==1
                  ylim([0.05,0.15]);
                  ylabel('Confidence Interverl Widths'); 
                  legend('IID FIM','IID Sobol','OU FIM','OU Sobol');
             end
             if j==2
                ylim([3,10]);
             end
             
              if j==3
                ylim([3.2,8]);
             end
        %         ylabel('Practical Identifiable Rate over 100 runs');
        %         legend('IID FIM','OU FIM','OU Sobol','Location','Northeast');

             title(sprintf('Parameter %s', prm_name{j}));   
        end
end
 





