function [P,Rate]=FP_algorithm(H,transmit_SNR,tolerance)

[Nt,K]=size(H);
noise=ones(K,1);

Pt=db2pow(transmit_SNR);

P_p=(0.4*Pt)/K; 

P=H./vecnorm(H)*sqrt(P_p);

flag=1;
obj_past=0;
count=0;

while (flag)


    T_k=sum(square_abs(H'*P),2)+noise;
    alpha_k=T_k./(T_k-square_abs(diag(H'*P)))-1;
    beta_k=sqrt(1+alpha_k).*diag(H'*P)./T_k;

    mu_min=0; %test wether mu=0
    warning('off','MATLAB:nearlySingularMatrix')
    P=((abs(beta_k').*H)*(abs(beta_k').*H)'+mu_min*eye(size(H,1)))\(transpose(sqrt(1+alpha_k).*beta_k).*H);
    warning('on','MATLAB:nearlySingularMatrix')
    if sum(sum(pow_abs(P,2)))-Pt<=0
        continue
    end

    % bisection search to find optimal dual variable mu
    muleft=0;
    muright=1;

    while 1

        P=((abs(beta_k').*H)*(abs(beta_k').*H)'+eye(Nt)*muright)\(transpose(sqrt(1+alpha_k).*beta_k).*H);
        
        if sum(sum(pow_abs(P,2)))-Pt<=0
            break
        end

        muright=muright*10;
    end

    while 1
        mu=(muleft+muright)/2;

        P=((abs(beta_k').*H)*(abs(beta_k').*H)'+eye(Nt)*mu)\(transpose(sqrt(1+alpha_k).*beta_k).*H);

       distance=sum(sum(pow_abs(P,2)))-Pt;

        if abs(distance)<=tolerance
            break
        end

        if distance>0
            muleft=mu;
        else
            muright=mu;
        end

    end

    T_ktest=sum(square_abs(H'*P),2)+noise;
    rate_p=log2(T_ktest./(T_ktest-square_abs(diag(H'*P))));
    obj=sum(rate_p);

    if abs(obj-obj_past)<=tolerance
        flag=0;
    else
        obj_past=obj;
        count=count+1;
    end
    if count>=800
        break;
    end
end

Rate=obj;
end




