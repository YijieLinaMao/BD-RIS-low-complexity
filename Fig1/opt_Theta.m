function out=opt_Theta(Q,lambda,gamma,A,a)
Nr=length(lambda);
R=length(lambda(lambda>1e-4));
x=zeros(Nr,Nr);
x=x+Q(:,R)*Q(:,R)';
for n=Nr-R+1:Nr
x=x+Q(:,n)*Q(:,n)'/(gamma-lambda(n));
end
 
out=x*A'*a;
end


