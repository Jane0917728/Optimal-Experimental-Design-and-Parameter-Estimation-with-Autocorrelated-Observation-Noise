   function ff=obj_con(r,K,phi,N0,x) %for both IID and OU
   
 
% r=0.2;K=50;sigma=1.5;phi=0.02;N0=1;
% x=linspace (0.1,80,6)';
%     N0=x(1);
%     x=t_best{3};
    t=x(1:end); % row vector
%     t=x10{1}';
    n_s=length(t);
    num_param=3;
   SIG_OU=zeros(n_s,n_s); % The covariance matrix of the Brownian bridge 
% The covariance matrix  of OU num_time*num_time 
    for s = 1:n_s
        for tt = 1:n_s
            SIG_OU(s, tt) = exp(-phi*abs(t(tt)-t(s)))/(2*phi);   
%            SIG_OU(s, tt) = Sigma^2*(exp(-phi*abs(t(tt)-t(s)))-exp(-phi*(t(tt)+t(s))))/(2*phi);    
        end   
    end 
   SIG=SIG_OU; %OU   
%   SIG=eye(n_s,n_s)*Sigma^2;%IID
   SIG_pinv=pinv(SIG);


%     mu_OU=x0*exp(-phi*t)+mu*(1-exp(-phi*t));
%         SIG_pinv=eye(n_s,n_s)/(Sigma^2);
  
pNpK=@(r,K,t,N0) N0^2*(1-exp(-r*t))./((K-N0)*exp(-r*t)+N0).^2;
pNpr=@(r,K,t,N0) K*N0*(K-N0)*exp(-r*t).*t./((K-N0)*exp(-r*t)+N0).^2;
pNpN0=@(r,K,t,N0) K^2*exp(-r*t)./((K-N0)*exp(-r*t)+N0).^2;
d_param=zeros(num_param,n_s);
param=[r;K;N0];
d_param(1,:)=pNpr(r,K,t,N0)*param(1);
d_param(2,:)=pNpK(r,K,t,N0)*param(2);
d_param(3,:)=pNpN0(r,K,t,N0)*param(3);
F=d_param*SIG_pinv*d_param';
 
%     ff=-det(FIM);
   ff=-log(det(F));
%         ff=-log(min(eig(FIM)));
  end