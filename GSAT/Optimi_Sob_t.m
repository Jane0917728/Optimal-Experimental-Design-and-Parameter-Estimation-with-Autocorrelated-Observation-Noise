%Main program for optimization of the number of the observed data time location or t using Sobol index

%Parameter definitions
num_time=30;
r=0.2; K=50;
lb_param=[0.17,44];
ub_param=[0.23,56];

num_points=8;
t_initial=0.1;
tf=60;N0=1;
fN=@(r,K,t,N0)  K*N0./((K-N0).*exp(-r.*t)+N0);
t=linspace(t_initial,tf,num_time)';
param_list=['r','K'];
rho=0.98;
phi=0.98;
num_MC=20000;% number of Monte Carlo

%  measurements=f(r,K,t,N0)+epsilon;
% SobrK=zeros(num_time,4);
% for i=1:num_time
%     SobrK(i,:)=Get_Sobol(param_list, lb_param,ub_param,N0,t(i),num_MC);
% end

S=Get_Sob_analyzed(r,K,lb_param,ub_param,t,N0);
dtheta_first=S;
% dtheta_first=[SobrK(:,1),SobrK(:,3)];
% dtheta_first(dtheta_first<0)=0; 
% dtheta_total=[SobrK(:,2),SobrK(:,4)];
% dtheta_total(dtheta_total<0)=0; 
 %compute the information for all num_time points according to the formular 
 % The element of  the FIM for each measurement point (t=1,...,num_time)
 Fim_correct=dtheta_first;
 for s=1:num_time
     for k=1:s-1
         for j=1:2
             Fim_correct(s,j)=Fim_correct(s,j)+(-1)^k*(rho+phi)*phi^(k-1)*dtheta_first(s,j);
         end
     end
 end
 
 Fim_global_ARMA=Fim_correct;
% Fim=zeros(num_time,2);
% for i=1:2 
%    Fim(1,i)=dtheta_first(1,i);
% end
% for tt=2:num_time  
%  for j=1:2
%    for i=1:tt-1           
%         Fim(tt,j)=Fim(tt,j)+(-phi)^(tt-i)*dtheta_first(i,j)-rho*(-phi)^(tt-i-1)*dtheta_first(i,j);
%    end
%    Fim(tt,j)=Fim(tt,j)+ dtheta_first(tt,j);
%  end
% end            

Fim_r_full=Fim_correct(:,1);
Fim_K_full=Fim_correct(:,2);

J_full=[Fim_r_full,Fim_K_full];
FIM_full=J_full'*J_full;

  eigFIM_full=eig(FIM_full);
  Sloppy_full=max(eigFIM_full)/min(eigFIM_full);
  Det_full=det(FIM_full);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lower and upper bounds
lb=zeros(num_time,1)';
ub=ones(num_time,1)';
% intcon=1:num_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 f1=@(x) obj_D(x,Fim_correct);
 f2=@(x) obj_S(x,Fim_correct);

%constraint
g1=@(x) sum(x(1:num_time))-num_points;
g2=@(x) -sum(x(1:num_time))+num_points;
%encoding
enc_b = zeros(1,num_time);
% enc_c = [1 1];
% enc = [enc_b enc_c];
% enc=[enc_b,1,1];
enc=enc_b+2;

%platemo('save',1,'algorithm',@SparseEA,'objFcn',{f},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub)
%   [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%  [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%     [Dec,Obj,Con] = platemo('algorithm',@MCCMO,'objFcn',{f1,f2},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
    [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1,f2},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);
%   [Dec,Obj,Con] = platemo('algorithm',@CMEGL,'objFcn',{f1},'encoding',enc,'conFcn',{g1,g2},'lower',lb,'upper',ub);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************************************************%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% find the dicision parameters from x

% tf=x(end)
% x=Dec(end,:);
% [obj_smalleast, obj_index] = min(Obj(:,2));
% x=Dec(obj_index,:);
% obj11=obj_smalleast

% [obj_smalleast, obj_index] = min(0.5*Obj(:,1)/max(Obj(:,1))+0.5*Obj(:,2)/max(Obj(:,2)));
 [obj_smalleast, obj_index] = min(Obj(:,1));
x=Dec(obj_index,:);
obj_com=obj_smalleast
obj11=obj_D(x,Fim_correct)
obj21=obj_S(x,Fim_correct)
% obj21=obj2_Sob(x,Fim)
 
t_index=find(x==1); 
t_opt=t(t_index); 
t1=linspace(t_initial,tf,num_time*2)';
 
N_selection=fN(r,K,t_opt,N0);
N=fN(r,K,t1,N0);
figure;
plot(t1, N, '-',t_opt,N_selection,'o','LineWidth', 2);  
xlabel('time');
ylabel('Density');
legend('Model curve', 'Selected points');
title('ARMA(1,1) noise with \rho=0.98, \phi=0.98 using the global Sobol index')


 