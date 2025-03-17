function N=Nfunction3(x,t)
r=x(1);
K=x(2);
N0=x(3);
% N=N0.*K.*exp(r.*t)./(K+N0.*(exp(r.*t)-1));
N=N0.*K./((K-N0)*exp(-r.*t)+N0);
end