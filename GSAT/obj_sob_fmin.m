  function ff= obj_sob_fmin(r,K, lb_param,ub_param,n,x)
% r=0.1; K=1000; w=0.8; tf=100; nmax=50; x_initial=[30;5;zeros(nmax-1,1)];
% x=x_initial;
% Maximize the objective function, so one doesn't need to find the inverse
% of a matrix.
 n=8;
N0=x(1);
tf=x(end);
Delta_t=x(2:n+1); 
%Delta_t(end)=tf-sum(Delta_t(1:n-2)); 
%%%%%%%%%%%%%%%******************************
% w2=(1-w)/2;
% w3=(1-w)/2;
%%%%%%%%%%%%%%%******************************
%find n from the dicision variable vector x
t=zeros(n,1);
for i=1:n
    t(i)=sum(Delta_t(1:i));
end
S=Get_Sob_analyzed(r,K,lb_param,ub_param,t,N0);

J_fmin=S; 
CV=J_fmin'*J_fmin/n; %partial N/partial (Log(r))=

% CV=CV/n;
% FIM=CV;
% CV=CV/sigma^2;

%D-optimal
ff= -det(CV);


%% E-optimal 
% eigCV=eig(CV);
% ff=min(eigCV);
% 
% 

%S-optimal
% eigCV=eig(CV);
% ff=min(eigCV)/max(eigCV);
% % % 


%%%Trace
%  eigCV=eig(CV);
%  ff=sum(1./eigCV);
%  ff=log(ff);

% ff=-log(ff);
%  ff=-log(ff);
% %***************************************


  end