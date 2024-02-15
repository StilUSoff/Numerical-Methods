format long

 folder = ['matrix_res/'];
 fileID = fopen([folder 'matrix_size.txt'],'w');
for n=5:53
    G=numgrid('B',n);
	A=delsq(G);
    s=length(A);
     writematrix(A,[folder 'A_matrix_' int2str(n) '.txt'],'Delimiter','tab')
     nbytes = fprintf(fileID,' %4s',int2str(n));
     nbytes = fprintf(fileID,' %4s\n',int2str(s));
    disp(n)
end
 fclose(fileID);