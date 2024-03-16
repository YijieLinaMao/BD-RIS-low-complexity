function out= Projection_of_Proposed(G,H,E)
%PROJECTION_OF_PROPOSED calculate the F-norm of the proposed solutions with
% various architectures
%   Input parameters:
% G----direct link between BS and users.
 % H----link between RIS and users
 % E----link between BS and RIS
% output parameters:
% F-norm of the three architectures

for Type_BD_RIS=1:3
    switch Type_BD_RIS
        case 1 %fully connected
            Theta=symuni(H*G'*E');
            F=G+E'*Theta'*H;
            out(1)=sqrt(trace(F*F'));
        case 2 % single connected
            Theta=symuni(diag(diag(H*G'*E')));
            F=G+E'*Theta'*H;
            out(2)=sqrt(trace(F*F'));
        case 3 % group connected with N_g=4;
            Theta=group_symuni(H*G'*E',4);
            F=G+E'*Theta'*H;
            out(3)=sqrt(trace(F*F'));
    end
end
end

