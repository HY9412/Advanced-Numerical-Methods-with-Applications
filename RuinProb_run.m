clear all;  clc;
 
phi=1.5; 
mu_X=1;
q=discretePareto(phi,mu_X);

u=10;                          %inithial reserve
gammaPI=0.95;                  %Poisson intensity
[kappa,psi]=RuinESM(u,gammaPI,q)

