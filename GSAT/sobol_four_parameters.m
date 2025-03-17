%% Sobol主程序

%若用于复杂模型，kp为全样本空间，将kp代入模型计算，根据输入和输出构造函数Sobol_obj。
clear all; close all; clc;
n_p = 4;  % 待分析参数数目
PS=[]; % 参数样本空间（一）
comp_PS=[]; % 参数样本空间（二），用于Monte Carlo采样
n_base= 1000; % 参数样本数目
N = n_base*(n_p*2+1); % 模型/函数运行总次数

PS=[ceil((30000 + 20000.*rand(n_base,1))) ceil((30000 + 20000.*rand(n_base,1))) ceil((0.8+ 0.4.*rand(n_base,1))) ceil((4+ 5.*rand(n_base,1)))]; % 对参数随机抽样
comp_PS=[ceil((30000 + 20000.*rand(n_base,1))) ceil((30000 + 20000.*rand(n_base,1))) ceil((0.8+ 0.4.*rand(n_base,1))) ceil((4+ 5.*rand(n_base,1)))]; % 对参数随机抽样
output=[];
c_out_1=[];
c_out_t=[];

% 计算模型输出
t=0;
for i=1:n_base
   t=t+1;
   kp(t,:)=PS(i,:);
   output(i,:)=Sobol_obj(kp(t,:)); %代入目标函数/模型计算
   for j=1:n_p
      t=t+1;
      kp(t,:)=[comp_PS(i,1:j-1),PS(i,j),comp_PS(i,j+1:n_p)]; % 构造Sobol抽样
      c_out_1(i,:,j)=Sobol_obj(kp(t,:)); %代入目标函数/模型计算
    
      t=t+1;
      kp(t,:) = [PS(i,1:j-1),comp_PS(i,j),PS(i,j+1:n_p)]; % 构造Sobol抽样
      c_out_t(i,:,j)=Sobol_obj(kp(t,:)); %代入目标函数/模型计算
   end
end
% t=N here;
%蒙特卡洛积分
n_out = size(output,2); 
f0 = zeros(1,n_out); 
D=zeros(1,n_out);  

for i = 1:n_base
 f0 = f0+output(i,:)/n_base; % 模型积分
 D = D+output(i,:).^2/n_base; % 计算模型输出方差
end
 
D=D-f0.^2; %模型输出方差
Dj=ones(n_p,1)*D;  
Dtotj=zeros(n_p,1); 
for i = 1:n_base
 for j = 1:n_p
 Dj(j,:)=Dj(j,:)-(output(i,:)-c_out_1(i,:,j)).^2/(2*n_base); %计算偏方差
 Dtotj(j,:)=Dtotj(j,:)+(output(i,:)-c_out_t(i,:,j)).^2/(2*n_base); %计算参数j的总方差
 end
end
%计算敏感度
Sob_1 = Dj./(ones(n_p,1)*D); %first order effect 一阶敏感度
Sob_t = Dtotj./(ones(n_p,1)*D); % total effect 总敏感度
disp(Sob_1); %负值或总和大于1为截断误差、非线性模型等因素造成
disp(Sob_t); %负值或总和大于1为截断误差、非线性模型等因素造成
Sob=[Sob_1,Sob_t];
bar(Sob)
xlabel('Variable index')
ylabel('Sobol index')
legend('First-order sensitivity','Total sensitivity','Best')
set(gca,'XTickLabel',{'kr','kf','ε','Ω'});

%%
ans1=Sobol_obj(kp(8,:));
%% 样例函数，需单独建立一个文件

function out=Sobol_obj(kp)
% a function that works with the trial Sobol Method program
% the trial function is fx: x1+5*x2*x3^2;
kr=kp(1);
kf=kp(2);
e=kp(3);
omiga=kp(4);
out= (kf+kr)./(2.*omiga.^2)+(kf.^3+kr.^3)./(2.*omiga.^2.*e.*kr.*kf)-(((e.*kf.^2.*kr-kf.^3).^2+(e.*kr.^2.*kf-kr.^3).^2-2.*e.*(e.*kf.^3.*kr.^3-kf.^2.*kr.^4-kr.^2.*kf.^4)).^0.5)./(2.*omiga.^2.*e.*kf.*kr);
end
