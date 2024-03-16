function sum_rate = Alternative_optimization_Manifold(G,H,E,tolerance,SNR)
%ALTERNATIVE_OPTIMIZATION_MANIFOLD 此处显示有关此函数的摘要
%   此处显示详细说明
H_d=G;
H_r=H;
G=E;
Pt=db2pow(SNR);
[~,K]=size(G);
Nr=size(H_r,1);
noise=ones(K,1);
%initial Theta as indentity matrix
Theta=eye(Nr);

H=H_d+G'*Theta*H_r;
P_p=(0.4* Pt)/K;
P=H_d./vecnorm(H_d)*sqrt(P_p);
 
flag=1;
t_past=0;
count=0;

while(flag)

    T_k=sum(square_abs(H'*P),2)+noise;
    alpha_k=T_k./(T_k-square_abs(diag(H'*P)))-1;
    beta_k=sqrt(1+alpha_k).*diag(H'*P)./T_k;

    res_k=log2(1+alpha_k)-alpha_k; 

    [P,~]=FP_update_P(H,Pt,alpha_k,beta_k,tolerance,res_k);

    [Theta,t]=Manopt_update_Theta(H_d,H_r,G,P,alpha_k,beta_k,Theta,noise);

    H=H_d+G'*Theta'*H_r;

    if abs(t-t_past)<=tolerance
        flag=0;
    else
        t_past=t;
        count=count+1;
    end
    if count>=2000
        break;
    end
end

%project Theta onto the symmetric unitary matrix
Theta=symuni(Theta);
F=H_d+G'*Theta'*H_r;

%calculate the final sum_rate
[~,sum_rate]=FP_algorithm(F,SNR,tolerance);

end

