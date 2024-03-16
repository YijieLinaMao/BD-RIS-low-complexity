function [G,H,E]=channel(seed,k,K,Nt,Nr)
 % Input parameters:
 % seed------control the random seed for the rng function
 % k-----sometimes we want to generate k channel realizations at one time.
 % K-----number of users
 % Nt----number of transmit antennas
 % Nr----number of BD-RIS elements 
 % Output parameters:
 % G----direct link between BS and users.
 % H----link between RIS and users
 % E----link between BS and RIS
 % We generate channels in this function. Specifically, the channels are
 % normalized such that the noise equal to 1.

    mu=inv(sqrt(10)*1e5); % normalization parameter

    alpha_d=3.5;
    alpha_r=2;
    alpha_g=2.2;  % path loss exponents

    T_0=-30; % reference pass loss at 1m

    rng(seed,"twister");
    IRS_location=50+1i*50; 

    user_location=150.+20*rand(K,1,k).*exp(1j*2*pi*rand(K,1,k));

    d_d=abs(user_location);
    d_r=abs(IRS_location);
    d_g=abs(user_location-IRS_location*ones(K,1,k)); %distance of each link


    L_d=sqrt(10^(T_0/10).*d_d.^(-alpha_d));
    L_r=sqrt(10^(T_0/10).*d_r.^(-alpha_r));
    L_g=sqrt(10^(T_0/10).*d_g.^(-alpha_g)); %large scale path loss


    
    G=1/sqrt(2)*(randn(Nt,K,k)+1i*randn(Nt,K,k));
    H=1/sqrt(2)*(randn(Nr,K,k)+1i*randn(Nr,K,k));
    E=1/sqrt(2)*(randn(Nr,Nt,k)+1i*randn(Nr,Nt,k));


    for i=1:k
       G(:,:,i)=G(:,:,i)*diag(L_d(:,:,i))./mu;
       H(:,:,i)=H(:,:,i)*diag(L_g(:,:,i))./sqrt(mu);
    end
    E=L_r*E./sqrt(mu);

end