function [alpha_Erl,T,Mean_PH]=erlang(k,lambda)
alpha_Erl=zeros(1,k);
alpha_Erl(1)=1;   %Initial distribution alpha.
T=zeros(k,k);
v=-lambda*ones(1,k);
u=lambda*ones(1,k-1);
T=diag(v)+diag(u,1);
Mean_PH=alpha_Erl*(-T)^(-1)*ones(k,1);

