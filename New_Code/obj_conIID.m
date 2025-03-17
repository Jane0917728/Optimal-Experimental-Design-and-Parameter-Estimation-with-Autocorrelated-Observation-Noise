   function ff=obj_conIID(r,K,N0,x) % The objective function with FIM for IID noise 
    t=x(1:end); % row vector
    n_s=length(t);
    num_param=3;    
  SIG_pinv= eye(n_s,n_s);


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
 
   ff=-log(det(F));
%         ff=-log(min(eig(FIM)));
  end