function lab


!polyfit research


x=-7:0.1:7;
y=2*x.^3-x.^2-60*x+1;
y=y+randn(size(y));

p1=polyfit(x,y,1);
p2=polyfit(x,y,2);
p3=polyfit(x,y,3);
p4=polyfit(x,y,4);

y1=polyval(p1,x);
y2=polyval(p2,x);
y3=polyval(p3,x);
y4=polyval(p4,x);

error=y-y4;

a=figure('position', [250 250 975 700]);
subplot(3,2,1)
plot(x,y,'-',...
    'Color','black',...
    'LineWidth',2)
ylim([-314 218]);
title(['polyval=1'], 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('y(x)', 'FontSize', 12);
grid on
hold on
scatter(x,y1,'.','m')
legend('actual', 'polyval');

subplot(3,2,2)
plot(x,y,'-',...
    'Color','black',...
    'LineWidth',2)
ylim([-314 218]);
title(['polyval=2'], 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('y(x)', 'FontSize', 12);
grid on
hold on
scatter(x,y2,'.','m')
legend('actual', 'polyval');

subplot(3,2,3)
plot(x,y,'-',...
    'Color','black',...
    'LineWidth',2)
ylim([-314 218]);
title(['polyval=3'], 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('y(x)', 'FontSize', 12);
grid on
hold on
scatter(x,y3,'.','m')
legend('actual', 'polyval')

subplot(3,2,4)
plot(x,y,'-',...
    'Color','black',...
    'LineWidth',2)
ylim([-314 218]);
title(['polyval=4'], 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('y(x)', 'FontSize', 12);
grid on
hold on
scatter(x,y4,'.','m')
legend('actual', 'polyval');


subplot(3,2,5)
scatter(x,error,'.','m')
title(['error(x)'], 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('error', 'FontSize', 12);
hold off

saveas(a,'bin/graphics/res.png')

!cftool research


y(11)=y(11)*4;
y(32)=y(32)*(-3);
y(58)=y(58)*7;
y(72)=y(72)*5;
y(93)=y(93)*(-9);
y(101)=y(101)*6;
y(124)=y(124)*(-11);
y(138)=y(138)*6;

A = x.';
B = y.';
save 'bin/X.dat' A -ascii;    
load (['bin/X.dat']);   
save 'bin/Y.dat' B -ascii; 
load (['bin/Y.dat']);   


!cftool spline
x1=-7:0.1:7;
y1=2*x.^3-x.^2-60*x+1;

A1 = x1.';
B1 = y1.';
save 'bin/X1.dat' A1 -ascii;    
load (['bin/X1.dat']);   
save 'bin/Y1.dat' B1 -ascii; 
load (['bin/Y1.dat']);


!cftool nonlinear
x2=-7:0.1:7;
y2=2*x.^3-x.^2-60*x+1;
y2=y2+randn(size(y));

A2 = x2.';
B2 = y2.';
save 'bin/X2.dat' A2 -ascii;
load (['bin/X2.dat']);   
save 'bin/Y2.dat' B2 -ascii; 
load (['bin/Y2.dat']);


cftool
end