function [Xout, ERR, z] = yakobi_matlab(A, b, eps)
D=diag(diag(A));
E=eye(size(A));
 
if isempty(find(diag(A)~=0))error('0 in D'),end 
D1= inv(D);
B=E-D1*A;
if(norm(B)>=1)error('Incorrect norm of B'),end
Xin=rand(length(b),1);
Xout=Xin;
for i=1:length(b)
    Xin(i)=1;
end
z=1;
while (z<100)
    Xout=B*Xin+D1*b;
    ERR(z)=norm(Xin-Xout);
    if(norm(Xin-Xout)<=eps)break,end
    z=z+1;
    Xin=Xout;
end
    
 return