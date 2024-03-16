function [Theta,out]=RS_update_Theta(H_d,H_r,G,P,alpha_k,beta_k,Theta_init,noise)
[Nr,K]=size(H_r);

[A_k,C_k]=deal(zeros(Nr,Nr));

e_k=G*P;
t_k=H_d'*P;

eta_k=log2(1+alpha_k)-alpha_k-square_abs(beta_k).*(sum(square_abs(t_k),2)+noise);
lambda_k=sqrt(1+alpha_k).*conj(beta_k).*diag(t_k);

for k=1:K
    j_user=1:K;

    z_k=conj(t_k(k,j_user)).*e_k*repmat(H_r(:,k)',K,1);

    A_k=A_k+sqrt(1+alpha_k(k))*conj(beta_k(k))*e_k(:,k)*H_r(:,k)'-square_abs(beta_k(k))*z_k;
    
    C_k=C_k+square_abs(beta_k(k))*H_r(:,k)*H_r(:,k)';
   
end
B_k=e_k*e_k';

% -2*real(trace(Theta_init*A_k))+real(trace(Theta_init*B_k*Theta_init'*C_k))

manifold=stiefelcomplexfactory(Nr, Nr);
problem.M=manifold;

problem.cost=@(Theta) -2*real(trace(Theta*A_k))+real(trace(Theta*B_k*Theta'*C_k));
problem.egrad=@(Theta) 2*C_k*Theta*B_k-2*A_k';

warning('off', 'manopt:getHessian:approx');
options = struct();
options.verbosity = 0;

[Theta, xcost] = conjugategradient(problem,Theta_init,options);

out=-xcost+sum(eta_k+2*real(lambda_k));

end

