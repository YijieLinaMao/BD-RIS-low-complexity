function [obj,newx,Theta] = quasi_Newton_update_Theta(H_d,H_r,G,P,noise,x)
 
Nr=size(H_r,1);
x_init=x;

object=@(x) 1e10*objective_function_in_quasi_Newton(x,H_d,H_r,G,P,noise);

options=optimoptions(@fminunc,'Display','off');
[xval,fval]=fminunc(object,x_init,options);
x=xval;
obj=-fval/1e10;
newx=x;
X=vectorx2matrixX(x,Nr);
Theta=(1i*X+50*eye(Nr))\(1i*X-50*eye(Nr));
end

