format long

folder = ['matrix_res/'];
fileID = fopen([folder 'obs.txt'],'w');
n=10;

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,1.09,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_1.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_1.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_1.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,1.45,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_2.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_2.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_2.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,1.9,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_3.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_3.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_3.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,5.5,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_4.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_4.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_4.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,10,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_5.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_5.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_5.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,19,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_6.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_6.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_6.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,46,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_7.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_7.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_7.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,91,10));
AD=Q*R*Q';
writematrix(AD,[folder 'A_matrix_8.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_8.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_8.txt'],'Delimiter','tab');


A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,1.9,10));
AD=Q*R*Q'; 
issymmetric(AD)
writematrix(AD,[folder 'A_matrix_unsim_1.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_unsim_1.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_unsim_1.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,19,10));
AD=Q*R*Q';
issymmetric(AD)
writematrix(AD,[folder 'A_matrix_unsim_2.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_unsim_2.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_unsim_2.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,1.9,10));
AD=Q*R*Q';
AD = triu(AD) + triu(AD,1)';
issymmetric(AD)
writematrix(AD,[folder 'A_matrix_sim_1.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_sim_1.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_sim_1.txt'],'Delimiter','tab');

A = rand(n);
[Q,R]=qr(A);
R=R-diag(diag(R))+diag(linspace(1,19,10));
AD=Q*R*Q';
AD = triu(AD) + triu(AD,1)';
issymmetric(AD)
writematrix(AD,[folder 'A_matrix_sim_2.txt'],'Delimiter','tab');
[V,D] = eig(AD);
fileID2=fopen([folder 'numbers_sim_2.txt'],'w');
for i=1:n
    nbyte = fprintf(fileID2,' %20.16f',D(i,i));
end
writematrix(V,[folder 'vect_sim_2.txt'],'Delimiter','tab');

nnnn=[0.01 0.05 0.1 0.5 1 2 5 10 0.1 2];
for i=1:10
    nbyte = fprintf(fileID,' %20.16f',nnnn(i));
end