function sum_rate = Alternative_optimization_quasi_Newton(G,H,E,tolerance,SNR)
%ALTERNATIVE_OPTIMIZATION_QUASI_NEWTON 
%   

H_d=G;
H_r=H;
G=E;
Nr=size(H_r,1);
x=ones(Nr*(Nr+1)/2,1);
Theta=vectorx2matrixX(x,Nr);
Pt=db2pow(SNR);
K=size(H_d,2);
noise=ones(K,1);
P_p=(Pt)/K; 
H=H_d+G'*Theta'*H_r;

P=H./vecnorm(H)*sqrt(P_p);

flag=1;
obj_past=0;
count=0;
maxcount=1000;
tic
while(flag)
    
    
    T_k=sum(square_abs(H'*P),2)+noise;
    alpha_k=T_k./(T_k-square_abs(diag(H'*P)))-1;
    beta_k=sqrt(1+alpha_k).*diag(H'*P)./T_k;

    res_k=log2(1+alpha_k)-alpha_k; 

    [P,~]=FP_update_P(H,Pt,alpha_k,beta_k,tolerance,res_k);
    
    [obj,x,Theta]=quasi_Newton_update_Theta(H_d,H_r,G,P,noise,x);
    H=H_d+G'*Theta'*H_r;
    if abs(obj-obj_past)<=tolerance
        flag=0;
    else
        obj_past=obj;
        count=count+1;
    end
    
    if count>=maxcount
        break;
    end
    
    
end

T_ktest=sum(square_abs(H'*P),2)+noise;
rate_p=log2(T_ktest./(T_ktest-square_abs(diag(H'*P))));

sum_rate=sum(rate_p);

end

