function sum_rate = classical_RZF(G,SNR)
%CLASSICAL_RZF 

K=size(G,2);
noise=ones(K,1);
F=G;

Pt=db2pow(SNR);
D=F/(F'*F+1/Pt*eye(K));
D=D./vecnorm(D);
P=diag(sqrt(Pt/K)*ones(K,1));
W=D*P;


T_k=sum(square_abs(F'*W),2)+noise;
rate_p=log2(T_k./(T_k-square_abs(diag(F'*W))));
sum_rate=sum(rate_p);

end

