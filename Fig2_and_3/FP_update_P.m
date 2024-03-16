function [P,out]=FP_update_P(H,Pt,alpha_k,beta_k,tolerance,res_k)

[Nt,K]=size(H);

mu_min=0;

warning('off','MATLAB:nearlySingularMatrix')
P=((abs(beta_k').*H)*(abs(beta_k').*H)'+mu_min*eye(size(H)))\(transpose(sqrt(1+alpha_k).*beta_k).*H);
warning('on','MATLAB:nearlySingularMatrix')

switch sum(sum(pow_abs(P,2)))-Pt<=0
    case 0

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

    case 1
end
noise=ones(K,1);
T_kcvx=sum(square_abs(H'*P),2)+noise;
E_k=-square_abs(beta_k).*T_kcvx+2*sqrt(1+alpha_k).*real(conj(beta_k).*(diag(H'*P)))+res_k;
out=sum(E_k);
end
