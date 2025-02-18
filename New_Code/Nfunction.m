function C=Nfunction(r,K,t,C0)
C=C0.*K.*exp(r.*t)./(K+C0.*(exp(r.*t)-1));
end