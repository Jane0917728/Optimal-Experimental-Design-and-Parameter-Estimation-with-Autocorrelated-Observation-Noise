 function ff=obj_autocorrelated(r,K,num_time,rho,phi,x)
% r=0.1; K=1000; w=0.8; tf=100; nmax=50; x_initial=[30;5;zeros(nmax-1,1)];
% x=x_initial;
% Maximize the objective function, so one doesn't need to find the inverse
% of a matrix.
N0=x(end-1);
tf=x(end);

variable=x(1:num_time)';% it is a 0,1 vector, where 1 denotes this point is selected as the new number.

t=linspace(0,tf,num_time)';


% f=@(r,K,t,N0) K*N0./((K-N0).*exp(-r*t)+N0);
% measurements=f(r,K,t,N0)+epsilon;

pNpK=@(r,K,t,N0) N0^2*(1-exp(-r*t))./((K-N0)*exp(-r*t)+N0).^2;

pNpr=@(r,K,t,N0) K*N0*(K-N0)*exp(-r*t).*t./((K-N0)*exp(-r*t)+N0).^2;


dr=pNpr(r,K,t,N0);
dK=pNpK(r,K,t,N0);
dtheta=[dr,dK];

 %compute the information for all num_time points according to the formular
 
 % The element of  the FIM for each measurement point (t=1,...,num_time)
 % suppose t>=1, when t<1, all variable=0;
 
      Fim=zeros(num_time,2);
      for j=1:2
%          Fim(1,i)=dtheta(1,i)*theta(i);
         Fim(1,j)=dtheta(1,j);
      end
      for tt=2:num_time
           for j=1:2
               for i=1:tt
                 if i<2
                    Fim(tt,j)=Fim(tt,j)+(-phi)^(tt-i)*dtheta(i,j);
                 else if i>=2
                         Fim(tt,j)=Fim(tt,j)+(-phi)^(tt-i)*dtheta(i,j)-rho*(-phi)^(tt-i)*dtheta(i-1,j);
                       end
                 end
                end
          end
%            for i=2:tt
%             for j=1:2
%                 Fim(tt,j)=Fim(tt,j)-rho*(-phi)^(tt-i)*dtheta(i-1,j);
%             end
%            end
      end
              Fim_r=Fim(:,1).*variable;
              Fim_K=Fim(:,2).*variable;
              J=[Fim_r*r,Fim_K*K];
              FIM=J'*J;
              FIM=FIM/sum(variable);
% FIM=CV;
% CV=CV/sigma^2;

%D-optimal
ff=det(FIM);


%% E-optimal 
% eigCV=eig(CV);
% ff=min(eigCV);
% 
% 

%%S-optimal
% eigCV=eig(CV);
% ff=min(eigCV)/max(eigCV);
% % % 


%%%Trace
%  eigCV=eig(CV);
%  ff=sum(1./eigCV);
%  ff=log(ff);


 ff=-log(ff);
% %*************



end