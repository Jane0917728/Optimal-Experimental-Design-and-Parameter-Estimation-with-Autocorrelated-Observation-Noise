%   x=Dec(end,:);
[obj_smalleast, obj_index] = min(0.5*Obj(:,1)/min(Obj(:,1))+0.5*Obj(:,2)/min(Obj(:,2)));
x=Dec(obj_index,:);
% [obj_smalleast, obj_index] = min(Obj(:,2));
% x=Dec(obj_index,:);
% obj11=obj_smalleast
 t_index=find(x==1); 


t1=linspace(t_initial,tf,num_time)';
 
t=t1(t_index);
N_selection=fN(r,K,t,N0);
N=fN(r,K,t1,N0);
figure;
plot(t1, N, '-',t,N_selection,'o','LineWidth', 2);  
xlabel('time');
ylabel('Density');
legend('Model curve', 'Selected points');
 
 
% title('IID noise');
title('ARMA(1,1) noise with \rho=0.98, \phi=0.98')