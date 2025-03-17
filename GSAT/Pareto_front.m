% Plot Pareto front
% N_selection=fN(r,K,t,N0);
% N=fN(r,K,t1,N0);
figure;
plot(Obj(:,1), Obj(:,2), 'o','LineWidth', 2);  
xlabel('D_opt');
ylabel('Sloppy');
% % legend('Data', 'Selected points');
title('IID noise');



% 假设x和y包含你的数据
 

% 将x和y转换为适合fit函数的格式
xData = Obj(:,1); % 确保是列向量
yData = Obj(:,2);

% 使用fit函数进行拟合，指定拟合类型，例如 'poly2' 代表二次多项式
% [fitResult, gof] = fit(xData, yData, 'poly1');
 [fitResult, gof] = fit(xData, yData,'linear');
%    [fitResult, gof] = fit(xData, yData,'power');

% 绘制拟合结果
figure;
plot(fitResult, xData, yData);
legend('data', 'fitting curves');
% xlabel('f_1');
% ylabel('f_2');
xlabel('D_{opt}');
ylabel('Sloppy');
title('Pareto Front');
 
