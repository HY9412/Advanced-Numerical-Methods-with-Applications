function q=discretePareto(phi,mu_X)
%gprnd(k,sigma,theta), tail index (shape) parameter k, scale parameter
%sigma, and threshold (location) parameter, theta.
N=10^4;
k_pareto=1/phi; mu=0; theta_pareto=mu; sigma=(phi-1)*mu_X/phi;
format long
q=gpcdf(1:N,k_pareto,sigma,theta_pareto)-gpcdf(0:N-1,k_pareto,sigma,theta_pareto);