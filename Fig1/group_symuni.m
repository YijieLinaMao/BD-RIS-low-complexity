function [X]=group_symuni(X,Ng)
Nr=size(X,1);
X=kron(eye(Nr/Ng),ones(Ng,Ng)).*X;

X=symuni(X);

end
