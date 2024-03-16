function out = Directly_Design(G,H,E,tolerance)
%Directly_Design calculate the F-norm of the proposed solutions with
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
           Theta=FC_DD(G,H,E,tolerance);
           F=G+E'*Theta'*H;
           out(1)=sqrt(trace(F*F'));
        case 2 % single connected
            Theta=SC_DD(G,H,E,tolerance);
            F=G+E'*Theta'*H;
            out(2)=sqrt(trace(F*F'));
        case 3 % group connected with N_g=4;
            Theta=GC_DD(G,H,E,tolerance);
            F=G+E'*Theta'*H;
            out(3)=sqrt(trace(F*F'));
    end
end

end

