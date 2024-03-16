function sum_rate = objective_function_in_quasi_Newton(x,H_d,H_r,G,P,noise)

 Nr=size(H_r,1);
 X=zeros(Nr);
 k=1;
 for i=1:Nr
     for j=i:Nr
         X(i,j)=x(k);
         k=k+1;
     end
     for j=1:i-1
         X(i,j)=X(j,i);
     end
 end
 Theta=(1i*X+50*eye(Nr))\(1i*X-50*eye(Nr));
 H=H_d+G'*Theta'*H_r;
 
 T_k=sum(square_abs(H'*P),2)+noise;
 alpha_k=square_abs(diag(H'*P))./(T_k-square_abs(diag(H'*P)));
 beta_k=sqrt(1+alpha_k).*diag(H'*P)./T_k;

 rate_k=log2(1+alpha_k)-alpha_k+2*sqrt(1+alpha_k).*real(conj(beta_k).*diag(H'*P))-square_abs(beta_k).*T_k;
 sum_rate=-sum(rate_k);
 
 end
 

