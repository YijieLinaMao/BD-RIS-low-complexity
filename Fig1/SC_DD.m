function Theta = SC_DD(G,H,E,tolerance)
%SC_DD 此处显示有关此函数的摘要
%   此处显示详细说明

H_RT=G';
H_RI=H';
H_IT=E;

[U,~,~]=svd(H_RT);

g=U(:,1)';
[~,~,V]=svd(H_RT);
w=V(:,1);

flag=1;
obj_past=0;
count=0;
maxcount=1000;

while (flag)


h_RT=g*H_RT*w;
h_ri=g*H_RI/norm(g*H_RI);
h_it=H_IT*w/norm(H_IT*w);



d=-transpose(angle(h_ri))-angle(h_it);

D=diag(exp(1j*d));

Theta=exp(1j*angle(h_RT))*D;



F=H_RT+H_RI*Theta*H_IT;


[U,~,V]=svd(F);

w=V(:,1);
g=U(:,1)';

obj=square_abs(g*F*w);

    if abs(obj-obj_past)<=tolerance
        flag=0;
    else
        obj_past=obj;
        count=count+1;
    end
    if count>=maxcount
        break;
    end

end

end

