folder = ['matrix_res/'];
n=10;
for k=1:6
    J = load([folder int2str(k) 'matrix.txt']);
    A=rand(n);
    o=0;
    for i=1:n
        for j=1:n
            o=o+1;
            A(i,j)=J(o);
        end
    end
    Xreal=rand(n,1);
    for i=1:n
        Xreal(i)=i;
    end
    b=rand(n,1);
    for i=1:n
        for j=1:n
    		b(i)=b(i)+A(i,j)*Xreal(j);
        end 
	end 
    b=A*Xreal;
    eps=0.0000000000000001;
    [Xout, ERR, z] = yakobi_matlab(A, b, eps)
end