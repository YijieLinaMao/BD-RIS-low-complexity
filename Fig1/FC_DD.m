function Theta = FC_DD(G,H,E,tolerance)
%FC_DD 此处显示有关此函数的摘要
%   此处显示详细说明
Nr=size(H,1);
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



A=0.5*(h_ri'*h_ri)+0.5*transpose(h_ri'*h_ri)-0.5*(h_it*h_it')-0.5*transpose(h_it*h_it');

[U,S]=eig(A);

delta=sort(diag(S),'descend');
U=fliplr(U);

T=eye(Nr);
T(1,1)=sqrt(-delta(Nr-1)/(delta(1)-delta(Nr-1)));
T(Nr-1,1)=sqrt(delta(1)/(delta(1)-delta(Nr-1)));

T(2,2)=sqrt(-delta(Nr)/(delta(2)-delta(Nr)));
T(Nr,2)=sqrt(delta(2)/(delta(2)-delta(Nr)));

T(:,5:Nr)=T(:,3:Nr-2);

b1=zeros(Nr,1);
b2=zeros(Nr,1);
b1(1)=T(Nr-1,1);
b1(Nr-1)=-T(1,1);

b2(2)=T(Nr,2);
b2(Nr)=-T(2,2);

T(:,3)=sqrt(1/2)*b1+sqrt(1/2)*b2;
T(:,4)=sqrt(1/2)*b1-sqrt(1/2)*b2;


V=U*T;

d=-transpose(angle(h_ri*V))-angle(V'*h_it);

D=diag(exp(1j*d));

Theta=exp(1j*angle(h_RT))*V*D*V.';


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

