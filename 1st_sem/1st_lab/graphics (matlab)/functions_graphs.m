g=figure;
a=1;
b=2;
x=a:0.00001:b;
plot(x,f1(x));
xlabel('x', 'fontsize', 12);
ylabel('f(x)', 'fontsize', 12);
saveas(g, 'function_1','png')

g=figure;
a=1;
b=2;
x=a:0.00001:b;
plot(x,f2(x));
xlabel('x', 'fontsize', 12);
ylabel('f(x)', 'fontsize', 12);
saveas(g, 'function_2','png')
clear;

g=figure;
a=2;
b=3;
x=a:0.00001:b;
plot(x,f3(x));
xlabel('x', 'fontsize', 12);
ylabel('f(x)', 'fontsize', 12);
saveas(g, 'function_3','png')

g=figure;
a=19;
b=22;
x=a:0.00001:b;
plot(x,f3(x));
xlabel('x', 'fontsize', 12);
ylabel('f(x)', 'fontsize', 12);
saveas(g, 'function_3_norm','png')

g=figure;
a=19;
b=22;
x=a:1:b;
plot(x,f3(x));
xlabel('x', 'fontsize', 12);
ylabel('f(x)', 'fontsize', 12);
saveas(g, 'function_3_big','png')