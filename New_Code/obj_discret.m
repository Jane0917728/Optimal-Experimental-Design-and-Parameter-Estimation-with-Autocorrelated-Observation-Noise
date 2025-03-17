    function ff=obj_discret(SIG,Sti,x)

%  x=[  0     0     0     1     0     0     1     1     0     0     1     0     1     0     0     0     0     0     0     0     0 0 0     0     0     0     0     0     0     0     0     1];
 
     % it is a 0,1 vector, where 1 denotes this point is selected as the new number.
    index=find(x==1);
%     t=t1(t_index);
% t_opt=t(x==1);
% num_time=length(t_opt);
% %    num_param=size(Sti,1);
% SIG_OU=zeros(num_time,num_time); % The covariance matrix  of OU num_time*num_time 
% for s = 1:num_time
%     for tt = 1:num_time
% %         SIG_OU(s, tt) = Sigma^2*(exp(-phi*abs(t(tt)-t(s)))-exp(-phi*(t(tt)+t(s))))/(2*phi);      
% %         SIG_OU(s, tt) = Sigma^2*exp(-phi*abs(t(tt)-t(s)))/(2*phi);  
%           SIG_OU(s, tt) = exp(-phi*abs(t_opt(tt)-t_opt(s)))/(2*phi);  
%     end   
% end 
% SIG_opt=SIG_OU;
    SIG_opt=SIG(index,index);  
    SIG_opt_pinv=pinv(SIG_opt);
    S_opt=Sti(:,index); 
%     F=zeros(num_param,num_param);
    F=S_opt*SIG_opt_pinv*S_opt';
%     for i=1:num_param
%           for j=i:num_param
%                 F(i,j)=S_opt(i,:)*SIG_opt_pinv*S_opt(j,:)';
%                 F(j,i)=F(i,j);
%           end
%        end
%   
%    ff=det(F_full)/det(F); 
  ff=-log(det(F));
% ff=-det(F);
% ff=-log(min(eig(F)));
 end