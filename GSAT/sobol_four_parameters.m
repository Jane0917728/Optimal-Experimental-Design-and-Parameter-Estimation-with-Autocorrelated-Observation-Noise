%% Sobol������

%�����ڸ���ģ�ͣ�kpΪȫ�����ռ䣬��kp����ģ�ͼ��㣬���������������캯��Sobol_obj��
clear all; close all; clc;
n_p = 4;  % ������������Ŀ
PS=[]; % ���������ռ䣨һ��
comp_PS=[]; % ���������ռ䣨����������Monte Carlo����
n_base= 1000; % ����������Ŀ
N = n_base*(n_p*2+1); % ģ��/���������ܴ���

PS=[ceil((30000 + 20000.*rand(n_base,1))) ceil((30000 + 20000.*rand(n_base,1))) ceil((0.8+ 0.4.*rand(n_base,1))) ceil((4+ 5.*rand(n_base,1)))]; % �Բ����������
comp_PS=[ceil((30000 + 20000.*rand(n_base,1))) ceil((30000 + 20000.*rand(n_base,1))) ceil((0.8+ 0.4.*rand(n_base,1))) ceil((4+ 5.*rand(n_base,1)))]; % �Բ����������
output=[];
c_out_1=[];
c_out_t=[];

% ����ģ�����
t=0;
for i=1:n_base
   t=t+1;
   kp(t,:)=PS(i,:);
   output(i,:)=Sobol_obj(kp(t,:)); %����Ŀ�꺯��/ģ�ͼ���
   for j=1:n_p
      t=t+1;
      kp(t,:)=[comp_PS(i,1:j-1),PS(i,j),comp_PS(i,j+1:n_p)]; % ����Sobol����
      c_out_1(i,:,j)=Sobol_obj(kp(t,:)); %����Ŀ�꺯��/ģ�ͼ���
    
      t=t+1;
      kp(t,:) = [PS(i,1:j-1),comp_PS(i,j),PS(i,j+1:n_p)]; % ����Sobol����
      c_out_t(i,:,j)=Sobol_obj(kp(t,:)); %����Ŀ�꺯��/ģ�ͼ���
   end
end
% t=N here;
%���ؿ������
n_out = size(output,2); 
f0 = zeros(1,n_out); 
D=zeros(1,n_out);  

for i = 1:n_base
 f0 = f0+output(i,:)/n_base; % ģ�ͻ���
 D = D+output(i,:).^2/n_base; % ����ģ���������
end
 
D=D-f0.^2; %ģ���������
Dj=ones(n_p,1)*D;  
Dtotj=zeros(n_p,1); 
for i = 1:n_base
 for j = 1:n_p
 Dj(j,:)=Dj(j,:)-(output(i,:)-c_out_1(i,:,j)).^2/(2*n_base); %����ƫ����
 Dtotj(j,:)=Dtotj(j,:)+(output(i,:)-c_out_t(i,:,j)).^2/(2*n_base); %�������j���ܷ���
 end
end
%�������ж�
Sob_1 = Dj./(ones(n_p,1)*D); %first order effect һ�����ж�
Sob_t = Dtotj./(ones(n_p,1)*D); % total effect �����ж�
disp(Sob_1); %��ֵ���ܺʹ���1Ϊ�ض���������ģ�͵��������
disp(Sob_t); %��ֵ���ܺʹ���1Ϊ�ض���������ģ�͵��������
Sob=[Sob_1,Sob_t];
bar(Sob)
xlabel('Variable index')
ylabel('Sobol index')
legend('First-order sensitivity','Total sensitivity','Best')
set(gca,'XTickLabel',{'kr','kf','��','��'});

%%
ans1=Sobol_obj(kp(8,:));
%% �����������赥������һ���ļ�

function out=Sobol_obj(kp)
% a function that works with the trial Sobol Method program
% the trial function is fx: x1+5*x2*x3^2;
kr=kp(1);
kf=kp(2);
e=kp(3);
omiga=kp(4);
out= (kf+kr)./(2.*omiga.^2)+(kf.^3+kr.^3)./(2.*omiga.^2.*e.*kr.*kf)-(((e.*kf.^2.*kr-kf.^3).^2+(e.*kr.^2.*kf-kr.^3).^2-2.*e.*(e.*kf.^3.*kr.^3-kf.^2.*kr.^4-kr.^2.*kf.^4)).^0.5)./(2.*omiga.^2.*e.*kf.*kr);
end
