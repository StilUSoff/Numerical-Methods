prompt = "what function? ";
num = input(prompt)
l = load('iterflcut_sec.txt')
if num==2
    f(num)=l(9,1);
    disp(f(num));
else
    f(num)=l(8,1);
    disp(f(num));
end

g=figure;
l = load('iterflcut_bis.txt');
x = l(:,2);
y1 = f(num)-l(:,1);
scatter(x,log10(y1), "filled","markerfacecolor","#7decf0");
title('iterations', 'fontsize', 12);
xlabel('iter', 'fontsize', 12);
ylabel('x-root(f)', 'fontsize', 12);

l = load('iterflcut_sec.txt');
x = l(:,2); 
y2 = f(num)-l(:,1);
hold on; 
scatter(x,log10(y2), "filled", "markerfacecolor","#d16608");
legend('bis', 'sec')
hold off;
saveas(g,'iter.png')


g=figure;
l = load('roots_bis.txt');
y1 = abs(f(num)-l(:,1));
x = l(:,3); 
scatter(log10(x),log10(y1), "filled","markerfacecolor","#7decf0");
title('epsilon', 'fontsize', 12);
xlabel('log10(eps)', 'fontsize', 12);
ylabel('log10(x-root(f))', 'fontsize', 12);

l = load('roots_sec.txt');
y2 = abs(f(num)-l(:,1));
x = l(:,3); 
hold on;
scatter(log10(x),log10(y2), "filled",  "markerfacecolor","#d16608");

hold on;
plot(log10(x),log10(x));
legend('bis', 'sec', 'bisector')
hold off;
saveas(g, 'eps.png')


g=figure;
l = load('roots_bis.txt');
x = l(:,3); 
y1 = l(:,2);
scatter(log10(x),y1, "filled","markerfacecolor","#7decf0");
title('number of iterations on epsilon', 'fontsize', 12);
xlabel('log10(eps)', 'fontsize', 12);
ylabel('number of iterations', 'fontsize', 12);

l = load('roots_sec.txt');
x = l(:,3); 
y2 = l(:,2);
hold on;
scatter(log10(x),y2, "filled","markerfacecolor","#d16608");
legend('bis', 'sec')
hold off;
saveas(g, 'iter_on_eps.png')


if num==2
g=figure;
l = load('randfluct_bis.txt');
x = l(:,1); 
y1 = l(:,2);
scatter(x,y1*100, "filled","markerfacecolor","#7decf0");
title('fluctuation', 'fontsize', 12);
xlabel('fluct', 'fontsize', 12);
ylabel('% of accuracy', 'fontsize', 12);

l = load('randfluct_sec.txt');
x = l(:,1); 
y2 = l(:,2);
hold on;
scatter(x,y2*100, "filled","markerfacecolor","#d16608");
legend('bis', 'sec')
hold off;
saveas(g, 'fluct.png');
end