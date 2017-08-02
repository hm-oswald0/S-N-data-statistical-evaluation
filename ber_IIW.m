function [N,d_sigma,k,s_log_N] = ber_IIW(N_A,L_a)

log_N = log10(N_A);
log_d_sigma = log10(L_a);

log_N_reg = [ones(size(log_N)) log_N];
beta = log_N_reg\log_d_sigma;

n = length(log_N);

k = -1/beta(2);
logC = -beta(1)/beta(2);

x_i = log_N + k*log_d_sigma;

x_m = sum(x_i)/n;

Stdv = sqrt(sum((x_m-x_i).^2)/(n-1));

k_IIW = 1.645*(1+1/sqrt(n));

x_k = x_m-k_IIW*Stdv;
x_k2 = x_m+k_IIW*Stdv;

s_log_N = Stdv;
N = [10^4;5*10^6];
d_sigma =[10^((logC-log10(N(1))-k_IIW*Stdv)./k),10^((logC-log10(N(1)))./k),10^((logC-log10(N(1))+k_IIW*Stdv)./k);...
          10^((logC-log10(N(2))-k_IIW*Stdv)./k),10^((logC-log10(N(2)))./k),10^((logC-log10(N(2))+k_IIW*Stdv)./k)];

end