function OUrandom=Generate_exact_OU(phi,Sigma,t_opt)


% Parameters for the OU process
% phi =0.01;       % Rate of mean reversion
% mu = 0;          % Long-term mean
% Sigma = 1.5;     % Volatility
% x_0= 0;          % Initial value of the process
% % 
% % % Time points to evaluate
%  t_opt = linspace(0,80,20);
% t_opt=[t_initial,t_opt];
num=length(t_opt);
% Allocate array for results
% X=zeros(1,num);
Z=randn(1,num);
% Number of steps
%x_0 follows N(0,\Sigma^2/(2\phi);
X0=Z(1)*sqrt(Sigma^2/(2*phi));
W=zeros(1,num); 
for i=1:num-1     
    W(i+1)=W(i)+sqrt(exp(2*phi*t_opt(i+1))-exp(2*phi*t_opt(i)))*Z(i+1);
end
ex=exp(-phi*t_opt);
X=X0*ex+Sigma*ex.*W/sqrt(2*phi);    
%    plot(t_opt, X, 'vg','MarkerSize',6,'LineWidth', 1.2);
X1=zeros(size(X));
 X1(1)=X0; 
   %fourth method
   for i=1:num-1
       dt=t_opt(i+1)-t_opt(i);
   mean_X = X1(i) * exp(-phi * dt) ;
    variance_X = (Sigma^2 / (2 * phi)) * (1 - exp(-2 * phi * dt));
    
    % Generate random value from normal distribution
    X1(i+1) = mean_X + sqrt(variance_X) * Z(i+1);
   end
OUrandom=X1(1:end);
% figure;
% plot(t_opt, X, '+m-',t_opt, X1, 'dg','MarkerSize',6,'LineWidth', 1.2);
 end


 
