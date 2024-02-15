i_n1=200
xx1=zeros(i_n1,1);
xx2=zeros(i_n1,1);

for n=10:10:2000
%1)
    i_n=n/10
    A=eye(n);
    for i=2:n
        for j=1:n-1
            if i>j
                A(i,j)=-1;
            elseif i<j
                A(i,j)=0;
            elseif i==j
                A(i,j)=1;
            end
        end
    end
    for j=1:n-1
        A(j,n)=1;
    end
%2)
    xt=rand(n,1);
    b=A*xt;
%3)
    x1=A\b;
    xx1(i_n)=norm(x1-xt)/norm(xt);
%4)
    [Q,R,p] = qr(A,"econ","vector");   
    x2(p,:) = R\(Q\b);
    xx2(i_n)=(norm(x2-xt)/norm(xt));
end

color1=['#ffafaf'];
color2=['#afffff'];
n=10:10:2000;
clear x1
clear x2

l=figure;
line(n,xx1,...
    'LineWidth',4,...
    'Color',color1);
hold on;
title('error(n)', 'fontsize', 12);
legend('Gauss');
xlabel('n', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(l, ['Gauss.png']);

g=figure;
line(n,xx2,...
    'LineWidth',4,...
    'Color',color2);
hold on;
title('error(n)', 'fontsize', 12);
legend('QR');
xlabel('n', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(g, ['QR.png']);