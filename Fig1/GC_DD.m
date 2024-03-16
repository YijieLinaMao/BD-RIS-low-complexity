function Theta_all = GC_DD(G,H,E,tolerance)
%GC_DD 此处显示有关此函数的摘要
%   此处显示详细说明

Nr=size(H,1);
Ng=4;
group=Nr/Ng;

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


h_RT_all=g*H_RT*w;
h_ri_all=g*H_RI/norm(g*H_RI);
h_it_all=H_IT*w/norm(H_IT*w);

for gr=1:group
    h_RT=h_RT_all;
    h_ri=h_ri_all(1+Ng*(gr-1):Ng*gr);
    h_it=h_it_all(1+Ng*(gr-1):Ng*gr);

A=0.5*(h_ri'*h_ri)+0.5*transpose(h_ri'*h_ri)-0.5*(h_it*h_it')-0.5*transpose(h_it*h_it');

[U,S]=eig(A);

delta=sort(diag(S),'descend');
U=fliplr(U);

T=eye(Ng);
T(1,1)=sqrt(-delta(Ng-1)/(delta(1)-delta(Ng-1)));
T(Ng-1,1)=sqrt(delta(1)/(delta(1)-delta(Ng-1)));

T(2,2)=sqrt(-delta(Ng)/(delta(2)-delta(Ng)));
T(Ng,2)=sqrt(delta(2)/(delta(2)-delta(Ng)));

T(:,5:Ng)=T(:,3:Ng-2);

b1=zeros(Ng,1);
b2=zeros(Ng,1);
b1(1)=T(Ng-1,1);
b1(Ng-1)=-T(1,1);

b2(2)=T(Ng,2);
b2(Ng)=-T(2,2);

T(:,3)=sqrt(1/2)*b1+sqrt(1/2)*b2;
T(:,4)=sqrt(1/2)*b1-sqrt(1/2)*b2;


V=U*T;

d=-transpose(angle(h_ri*V))-angle(V'*h_it);

D=diag(exp(1j*d));

Theta=exp(1j*angle(h_RT))*V*D*V.';

Theta_all(1+Ng*(gr-1):Ng*gr,1+Ng*(gr-1):Ng*gr)=Theta;
end

F=H_RT+H_RI*Theta_all*H_IT;


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

