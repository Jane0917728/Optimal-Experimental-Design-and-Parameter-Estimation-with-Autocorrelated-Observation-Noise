

num_timenew=25;
 
  SIG_new=SIG(1:num_time,1:num_time);
%   SIG_pinv_new=pinv(SIG_new);
SIG_pinv_new=pinv(SIG_new);
 Sti_new=ones(3,num_time);
 p1=[0.8,0.1,0.9];
 p2=[0.1,0.9,0.05];
 num1=floor(num_time/2);
 num2=num_time-num1;
 for i=1:3
 Sti_new(i,1:num1)=Sti_new(i,1:num1)*p1(i);
 Sti_new(i,num1+1:end)=Sti_new(i,num1+1:end)*p2(i);
 end
% F_new=zeros(num_param,num_param);
F_new=Sti_new*SIG_pinv_new*Sti_new';
%   for i=1:num_param
%         for j=i:num_param
%             F_new(i,j)=Sti_new(i,:)*SIG_pinv_new*Sti_new(j,:)';
%             F_new(j,i)=F_new(i,j);
%         end
%   end
  detF_new=det(F_new);
  
    variables=x;% it is a 0,1 vector, where 1 denotes this point is selected as the new number.
    index=find(variables==1);
%     t=t1(t_index);

   num_param=size(Sti,1)-1;
    SIG_opt=SIG_new(index,index);    
    SIG_opt_pinv=pinv(SIG_opt);
    S_opt_new=Sti_new(:,index); 
    F_1=S_opt_new*SIG_opt_pinv*S_opt_new';
     
  
%    ff=det(F_full)/det(F); 
det_Fopt=det(F_1);