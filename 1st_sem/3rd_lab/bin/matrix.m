i=[1;10;100;1000;10000;100000];
!folder = ['bin/matrix_res/'];
folder = ['matrix_res/'];
for j=1:6
    A = rand(10);
    [u, d, v] = svd(A);
    d = eye(10);
    d(1,1) = i(j); 
    M = u*d*v
    M=M';
    fileID = fopen([folder int2str(j) 'matrix.txt'],'w');
    nbytes = fprintf(fileID,' %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f\n',M);
    fclose(fileID);
    type([folder int2str(j) 'matrix.txt'])
end