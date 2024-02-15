function matrix(size)
format long

folder = ['matrix_res/'];
for j=1:6
A = rand(size);
[u, d, v] = svd(A);
d = eye(size);
d(1,1) = 100000*d(1,1);  
M100 = u*d*v'; 
cond(M100) ;
M100D=M100;
    for i=1:size
        M100(i,i);
        sum(M100(i,[1:i-1 i+1:size]));
        if j==1
            M100D(i,i)=M100D(i,i)+i*sum(M100(i,[1:i-1 i+1:size]));
        elseif j==2
            M100D(i,i)=M100D(i,i)+i*i*sum(M100(i,[1:i-1 i+1:size]));
        elseif j==3
            M100D(i,i)=M100D(i,i)+i*i*i*sum(M100(i,[1:i-1 i+1:size]));
        elseif j==4
            M100D(i,i)=M100D(i,i)+i*i*i*i*sum(M100(i,[1:i-1 i+1:size]));
        elseif j==5
            M100D(i,i)=M100D(i,i)+i*i*i*i*i*sum(M100(i,[1:i-1 i+1:size]));
        elseif j==6
            M100D(i,i)=M100D(i,i)+i*i*i*i*i*i*sum(M100(i,[1:i-1 i+1:size]));
        end
    end
    condi(j)=cond(M100D)
    for i=1:size
        M100D(i,i) > sum(M100D(i,[1:i-1 i+1:size]))
    end
    M100D=M100D';
    fileID = fopen([folder int2str(j) 'matrix.txt'],'w');
    nbytes = fprintf(fileID,' %25.16f\n',M100D);
    fclose(fileID);
end
    fileID = fopen([folder 'size.txt'],'w');
    nbytes = fprintf(fileID,' %4s\n',int2str(size));
    fclose(fileID);
    fileID = fopen([folder 'cond.txt'],'w');
    nbytes = fprintf(fileID,' %25.16f\n',condi);
    fclose(fileID);
end