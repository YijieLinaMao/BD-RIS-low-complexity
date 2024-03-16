function out = Projection_of_Optimal(G,H,E,tolerance)
%PROJECTION_OF_OPTIMAL calculate the F-norm of the optimal solution of the relaxed problem with
% various architectures
% Input parameters:
% G----direct link between BS and users.
% H----link between RIS and users
% E----link between BS and RIS
% output parameters:
% F-norm of the three architectures

Nr=size(H,1);
A=kron(E.',H');
a=vec(G');

[Q,D]=eig(A'*A);
lambda=diag(D);

gamma=max(lambda)+5;

theta=opt_Theta(Q,lambda,gamma,A,a);
if theta'*theta>Nr
    gamma=5*gamma;
end


count=0;
gamma_max=gamma;
gamma_min=max(lambda)+0.1;

while(count<1000)

gamma=(gamma_max+gamma_min)/2;

theta=opt_Theta(Q,lambda,gamma,A,a);

dis=theta'*theta-Nr;

if abs(dis)<=tolerance
    break
end

if dis>0
    gamma_min=gamma;
else
    gamma_max=gamma;
end
count=count+1;

end

Theta=reshape(theta,[Nr,Nr]);

Theta_fully=symuni(Theta);
F=G+E'*Theta_fully'*H;
out(1)=sqrt(trace(F*F'));

Theta_group=group_symuni(Theta,4);
F=G+E'*Theta_group'*H;
out(3)=sqrt(trace(F*F'));

Theta_single=group_symuni(Theta,1);
F=G+E'*Theta_single'*H;
out(2)=sqrt(trace(F*F'));

end

