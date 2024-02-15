%1
disp('researc 1')
v = [1 -2 3]
disp('Euclidean norm of vector')
n = norm(v)
disp('1-norm vector')
n = norm(v,1)

X = [2 0 1; -1 1 0; -3 3 0]
disp('2-norm of matrix')
n = norm(X)
disp('1-norm of matrix')
n = norm(X,1)
disp('Frobenius norm')
n = norm(X,"fro")



%2
disp('researc 2')
A=[1 2 3; 4 5 6; 7 8 9]
B=[1*10.^8 2*10.^8 3*10.^8; 4*10.^8 5*10.^8 6*10.^8; 7*10.^8 8*10.^8 9*10.^8]
det(A)
det(B)
A1=rand(3,1);
B1=rand(3,1);
X1=A\A1
X2=B\B1



%3
disp('researc 3')
H1 = hilb(5);
H2 = hilb(10);
H3 = hilb(15);
X1=[1;2;3;4;5];
X2=[1;2;3;4;5;6;7;8;9;10];
X3=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
B1=H1*X1;
B2=H2*X2;
B3=H3*X3;
X1err=H1\B1;
X2err=H2\B2;
X3err=H3\B3;
B1err=H1*X1err;
B2err=H2*X2err;
B3err=H3*X3err;

disp('cond(5x5)')
cond(H1)
disp('cond(10x10)')
cond(H2)
disp('cond(15x15)')
cond(H3)

error1=abs(X1-X1err);
error2=abs(X2-X2err);
error3=abs(X3-X3err);
nevazka1=abs(B1-B1err);
nevazka2=abs(B2-B2err);
nevazka3=abs(B3-B3err);

k=figure;
loglog(X1,error1,'-bo',...
    'LineWidth',1,...
    'Color',"#004A43",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#004A43");
hold on;
loglog(X2,error2,'-bo',...
    'LineWidth',1,...
    'Color',"#EB404E",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#EB404E");
hold on;
loglog(X3,error3,'-bo',...
    'LineWidth',1,...
    'Color',"#FFAC3E",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#FFAC3E");
hold on;
title('Error', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('X(i)-Xerr(i)', 'fontsize', 12);
xlim([1,15]);
hold off;
legend('5x5', '10x10', '15x15');
saveas(k, 'graphics\error.png');
close

j=figure;
plot(X1,nevazka1,'-bo',...
    'LineWidth',1,...
    'Color',"#004A43",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#004A43");
hold on;
plot(X2,nevazka2,'-bo',...
    'LineWidth',1,...
    'Color',"#EB404E",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#EB404E");
hold on;
plot(X3,nevazka3,'-bo',...
    'LineWidth',1,...
    'Color',"#FFAC3E",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#FFAC3E");
hold on;
title('Nevazka', 'fontsize', 12);
xlabel('i', 'fontsize', 12);
ylabel('B(i)-Berr(i)', 'fontsize', 12);
xlim([1,15]);
hold off;
legend('5x5', '10x10', '15x15');
saveas(j,'graphics\nevazka.png');
close



%4
disp('researc 4')
tEnd=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
k=1;
for i=5:5:100
    tStart_chol= tic;
    A = rand(i);
    [u, d, v] = svd(A);
    d = eye(i);
    d(1,1) = 1; 
    M = u*d*v';
    M = M*M';
    cond(M);
    issymmetric (M);
    [~,p] = chol(M);
    R = chol(M);
    b = rand(i,1);
    x = R\(R'\b);
    tEnd_chol(k,:) = toc(tStart_chol);
    tStart_LU= tic;
    x = M\b;
    [L,U,P] = lu(M);
    y = L\(P*b);
    x = U\y;
    tEnd_LU(k,:) = toc(tStart_LU);
    tStart_qr= tic;
    [Q,R,p] = qr(M,0);
    x(p,:) = R\(Q\b);
    tEnd_qr(k,:) = toc(tStart_qr);
    k=k+1;
end
i=5:5:100;
j=figure;
loglog(i,tEnd_chol,'-bo',...
    'LineWidth',1,...
    'Color',"#f36429",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#f36429");
hold on;
loglog(i,tEnd_LU,'-bo',...
    'LineWidth',1,...
    'Color',"#1baf4e",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#1baf4e");
hold on;
loglog(i,tEnd_qr,'-bo',...
    'LineWidth',1,...
    'Color',"#2a2a2a",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#2a2a2a");
hold on;
legend('chol', 'LU', 'qr');
title('time(NxN)', 'fontsize', 12);
xlabel('NxN', 'fontsize', 12);
ylabel('time', 'fontsize', 12);
hold off;
saveas(j,'graphics\time.png');
close



%5
disp('researc 5')
disp('If A-scalar:')
A=2
B=[2;4;6;8;10]
C=A\B
disp('A\B is equivalent to A.\B')

disp('If A-square N matrix and B-matrix with N rows:')
A=[1 2 3; 4 3 1; 6 2 3]
B=[1;1;1]
C=A\B
disp('A\B is a solution to the equation A*x = B')

disp('If A- matrix MxN where M~N and B-matrix with M rows:')
A=[1 2 3; 4 3 1; 6 2 3; 4 5 8;]
B=[1;1;1;1]
C=A\B
disp('A\B is a least-squares solution to the system of equations A*x= B')



%6
disp('researc 6')
A=[1 0 0; 4 3 0; 6 2 3]
B=[1;1;1]
opts.UT = true;  opts.LT = false; opts.POSDEF = false;
disp('opts.UT = true;  opts.LT = false; opts.POSDEF = false;')
X=linsolve(A,B,opts)

A=[1 2 3; 0 3 1; 0 0 3]
B=[1;1;1]
opts.UT = false;  opts.LT = true; opts.POSDEF = false;
disp('opts.UT = false;  opts.LT = true; opts.POSDEF = false;')
X=linsolve(A,B,opts)




%7
disp('researc 7')
rng default
A = sprand(3,3,.5);
A = A'*A.*10
B = sum(A,2)
x0 = [0;0;0]
X_res=A\B
disp('change maxit')
for k = 1:5
    disp('k=');
    disp(k);
    tol = 0.001;
    maxit = k
    [x,flag,relres] = pcg(A,B,tol,maxit,[],[],x0);
    X(:,k) = x;
    Flag_k=flag
    Error= relres
end
disp('solutions from k=1 to k=5:');
disp(X);
disp('change tol')
for k = 1:5
    disp('k=');
    disp(k);
    tol = 0.001/k
    maxit = 100;
    [x,flag,relres] = pcg(A,B,tol,maxit,[],[],x0);
    X(:,k) = x;
    Flag_k=flag
    Error= relres
end
disp('solutions from k=1 to k=5:');
disp(X);
disp('change x0')
for k = 1:5
    disp('k=');
    disp(k);
    tol = 0.001;
    maxit = 100;
    [x,flag,relres] = pcg(A,B,tol,maxit,[],[],x0);
    X(:,k) = x;
    Flag_k=flag
    Error= relres
    x0 = [k;k;k]
end
disp('solutions from k=1 to k=5:');
disp(X);




%8
disp('researc 8')
A=rand(5)
X=ones(5,1)
B=A*X
tol = 0.0000000000000001
maxit = 100
x0 = [0.5;0.5;0.5;0.5;0.5];
[x,flag,relres] = pcg(A,B,tol,maxit,[],[],x0);
Flag_k=flag
Error= relres
disp('solutions from pcg');
disp(x);




%9
disp('researc 9')

A=rand(10);
[u,d,v]=svd(A);
d=eye(10);
d(1,1)=10;
A=u*d*v;
A=A*A';
X=ones(10,1);
B=A*X;
tol=0.0000000000000001;
[x,flag,relres,iter] = pcg(A,B,tol,2);
disp('relres, iter, flag')
disp(relres)
disp(iter)
disp(flag)


z=0;
for k = 1:5:51
    z=z+1;
    A=rand(50);
    [u,d,v]=svd(A);
    d=eye(50);
    d(1,1)=k;
    A=u*d*v;
    A=A*A';
    for i=1:50
        B(i)=i
    end
    tol=0.0000000000000001;
    [x,flag,relres,iter] = pcg(A,B,tol);
    Iterations(z)=iter;
end
k=1:5:51;
v=figure;
plot(k,Iterations,'-bo',...
    'LineWidth',1,...
    'Color',"#004A43",...
    'MarkerSize',3,... 
    'MarkerFaceColor',"#004A43");
hold on;
title('iter(k)', 'fontsize', 12);
xlabel('k', 'fontsize', 12);
ylabel('iter', 'fontsize', 12);
hold off;
saveas(v,'graphics\iter(k).png');
close