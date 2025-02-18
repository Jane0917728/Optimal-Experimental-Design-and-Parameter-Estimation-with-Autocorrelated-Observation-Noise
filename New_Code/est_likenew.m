function L=est_likenew(param_list,SIG,measurements,t_opt)
%
 
r=param_list(1);
K=param_list(2);
N0=param_list(3);
% mu_opt=mu(t_index);
% SIG_opt=SIG(t_index,t_index);
SIG_pinv=pinv(SIG);
% measure_opt=measurements(t_index);
T=length(t_opt);
real=Nfunction(r,K,t_opt,N0);
err=measurements-real;
 
a1=det(SIG);
a2=log(det(SIG));
Likelihood=-T/2*log(2*pi)-1/2*log(det(SIG))-1/2*err*SIG_pinv*err';
 
L=-Likelihood;

end

