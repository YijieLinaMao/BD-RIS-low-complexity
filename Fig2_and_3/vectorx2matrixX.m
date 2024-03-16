function out=vectorx2matrixX(x,Nr)
 
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
out=X;
end
