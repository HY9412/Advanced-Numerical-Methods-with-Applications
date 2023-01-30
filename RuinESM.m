function [kappa,psi]=RuinESM(u,gammaPI,q)
    N = length(q);
    M = 700;                  %Truncation level of the infinit series
    k=4;lambda=4;             %Parameters of Erlang distribution
    [~,T,~]=erlang(k,lambda); %Erlang distribution as a special case of Phase-type distributions
    theta=max(-diag(T));    
    c0 = gammaPI / lambda; c1 = gammaPI / theta; 
    kappa = zeros(1,M+1);
    kappa(1) = (c0 * k) * sum((1:N).*q) ;
    ss = (lambda / theta) ./ (1:N);            
    counter=0;  sum3=zeros(1,M+1); CC=zeros(k,M); DD=zeros(k,M);
    
    %Calculate Ci and Di. Use parfor to speed up the calculations.
    tic
    parfor n=k:M            % n is the index of the kappa, for n<k we don't need to calculate the cdf of the Binomial distribution
            CC(n) = sum(q .* binocdf(k-2,n-1,ss));
            DD(n) = sum((1:N) .* q .* binocdf(k-1,n,ss));    
    end
    toc

    for n=1:M % n is the index of the kappa
        if n<k
            kappa(n+1) = (kappa(n)-1)*(c1+1)+1;
        elseif n==k
            sum11 = c1 * sum(kappa(1:k));  
            sum22 = c0 * k * DD(k) - n * c1 * CC(k);
            kappa(n+1) = sum11 + sum22;
        elseif n>k            
            sum1(n+1) = c1 * sum(kappa(n+1-k:n)); 
            sum2(n+1) = c1 * sum(kappa(1:n-k).*CC(1:n-k));            
            sum3(n+1) = c0 * k * DD(n) - n * c1 * CC(n);
            kappa(n+1) = sum1(n+1) + sum2(n+1) + sum3(n+1); 
        end          
%         fprintf('Just finished iteration #%d\n', counter);
        counter = counter + 1;
    end
    format long    
    figure(1)
    plot([0:M],kappa,'r*')
    figure(2)
    plot([0:M],sum1,'ko')
    hold on
    plot([0:M],sum2,'b.')
    plot([0:M],sum3,'r*')
    hold off
    figure(3)
    plot([0:M],sum1+sum2,'b.')
    hold on
    plot([0:M],kappa,'r*')
    hold off
    
    Psi = kappa.*poisspdf(0:M,theta*u); %The n-th term of the series is the product of kappa(n) with pdf of the Poisson distribution 
    psi=sum(Psi);
    
    
