tabl = zeros(14,10);
tabl2 = zeros(4,5);
tabl3 = zeros(14,2);
% а=1,2 б=3,4 в=5,6 г=7,8 д=9,10 е=11,12 ж=13,14


% !!! 1)

for i = 1:14
    if (i == 6)
        n = 200;
    elseif (i ==1) || (i == 5) || (i ==7) || (i == 11) || (i == 13)
        n = 100;
    else 
        n = 500;
    end
    A = rand(n);
    C = rand(n);
    [Q,R] = qr(A);
    if (i == 1) || (i ==2) 
        R = R-diag(diag(R))+diag(linspace(1,3*n,n));
    elseif (i == 3)
        R = R-diag(diag(R))+diag(linspace(1,3*n,n));
        for j = 1:n
            for m = 1:n
                if (R(j,m) == 0)
                    R(j,m) = 1e-16;
                end
            end
        end
    elseif (i == 4)
        R = R-diag(diag(R))+diag(linspace(1,3*n,n));
        for j = 1:n
            for m = 1:n
                if (R(j,m) == 0)
                    R(j,m) = 1e-100;
                end
            end
        end
    elseif (i == 5) || (i == 6)
        R = R-diag(diag(R))+diag(linspace(0,1,n));
    elseif (i == 7) || (i == 8)
        R = R-diag(diag(R))+diag(linspace(0,1,n));
        R(1,1) = 10;
    elseif (i == 9)
        R = R-diag(diag(R))+diag(linspace(0,1,n));
        R(1,1) = 10;
        for j = 1:n
            for m = 1:n
                if (R(j,m) == 0)
                    R(j,m) = 1e-16;
                end
            end
        end
    elseif (i == 10)
        R = R-diag(diag(R))+diag(linspace(0,1,n));
        R(1,1) = 10;
        for j = 1:n
            for m = 1:n
                if (R(j,m) == 0)
                    R(j,m) = 1e-100;
                end
            end
        end
    elseif (i == 11) || (i == 12)
        R = R-diag(diag(R))+diag(linspace(0,1,n));
        R(1,1) = -10;
    else
        R = R-diag(diag(R))+diag(linspace(0,3*n,n));
        R(1,1) = complex(1,2);
        R(2,2) = complex(1,-2);
        R(3,3) = complex(1,3);
        R(4,4) = complex(1,-3);
        R(5,5) = complex(1,4);
        R(6,6) = complex(1,-4);
    end 

    B = Q*R*Q';
    f = figure;
    spy(B, 'm')
    saveas(f, ['results/(' int2str(i) ')_matrix.png']);
    close


    % !!! 2)

    [V,D] = eig(B);
    disp(i)
    zero_check=norm(B*V-V*D)


    % !!! 3)

    tic;
    d = eig(B);
    Elapsed = toc;
    tabl(i,1) = Elapsed;

    tic;
    d = eig(B,C);
    Elapsed = toc;
    tabl(i,2) = Elapsed;

    tic;
    [V,D] = eig(B);
    Elapsed = toc;
    tabl(i,3) = Elapsed;
    eps = norm(B*V-V*D);
    tabl(i,4) = eps;

    tic;
    [V,D] = eig(B, 'nobalance');
    Elapsed = toc;
    tabl(i,5) = Elapsed;
    eps = norm(B*V-V*D);
    tabl(i,6) = eps;

    tic;
    [V,D] = eig(B,C);
    Elapsed = toc; 
    tabl(i,7) = Elapsed;
    eps = norm(B*V-V*D);
    tabl(i,8) = eps;

    tic;
    [V,D] = eig(B,C, 'chol');
    Elapsed = toc;
    tabl(i,9) = Elapsed;
    eps = norm(B*V-V*D);
    tabl(i,10) = eps;


    % !!! 4)

    if (i == 9) || (i == 10)
        F = hess(B);
        if (i == 9)
            l = 0;
        else
            l = 2;
        end
        for p = 1:2
            if (p == 1)
                F(n,1) = 1e-16;
            else
                F(n,1) = 1e-100;
            end
            tic;
            [V, D, flag] = eigs(F, 2, 'largestabs');
            Elapsed = toc;
            tabl2(p+l,1) = Elapsed;
            tic;
            [V, D, flag] = eigs(F, 2, 'smallestabs');
            Elapsed = toc; 
            tabl2(p+l,2) = Elapsed;
            tic;
            [V, D, flag] = eigs(F, 2, 'largestreal');
            Elapsed = toc;
            tabl2(p+l,3) = Elapsed;
            tic;
            [V, D, flag] = eigs(F, 2, 'smallestreal');
            Elapsed = toc; 
            tabl2(p+l,4) = Elapsed;
            tic;
            [V, D, flag] = eigs(F, 2, 'bothendsreal');
            Elapsed = toc; 
            tabl2(p+l,5) = Elapsed;
        end
    end


    % !!! 5)

    tic;
    H = hess(B);
    err = inf;
    k = 0;
    while err > 1e-2
        k = k+1;
        [Q1, R1] = qr(H);
        H = R1 * Q1;
        err = norm(tril(H, -1), 'fro');
        if (k>=5000)
            break
        end
    end
    Elapsed = toc;
    tabl3(i,1) = Elapsed;
    tabl3(i,2) = err;

end

n3=['3.1 time' '3.2 time' '3.3 time' '3.3 err' '3.4 time' '3.4 err' '3.5 time' '3.5 err' '3.6 time' '3.6 err'];
n4=['largestabs  ' 'smallestabs ' 'largestreal ' 'smallestreal' 'bothendsreal' ];
n5=['time' 'error'];

writematrix(tabl,['results/res_3.txt'],'Delimiter',' ');
writematrix(tabl2,['results/res_4.txt'],'Delimiter',' ');
writematrix(tabl3,['results/res_5.txt'],'Delimiter',' ');

