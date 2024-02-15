%folder = ['bin/researches_res/'];
folder = ['researches_res/'];
%folders = ['bin/graphics/'];
folders = ['graphics/'];
color=['#e6b3ff'];

J = load([folder 'error.txt']);
K = load([folder 'time.txt']);
D = load([folder 'fluct.txt']);
M = load([folder 'error_eps.txt']);
N = load([folder 'iter_eps.txt']);   
O = load([folder 'error_iter.txt']);
OO = load([folder 'error_p.txt']);


j=figure;
x = J(:,1);
y = J(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('error on matrix size', 'fontsize', 12);
xlabel('matrix size', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(j, [folders 'error.png']);
close
    
k=figure;
x = K(:,1);
y = K(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('time on matrix size', 'fontsize', 12);
xlabel('matrix size', 'fontsize', 12);
ylabel('time, c', 'fontsize', 12);
hold off;
saveas(k, [folders 'time.png']);
close

d=figure;
x = D(:,1)';
y = D(:,2)';
x0=1:1:3;
y0=[0 0 0];
x1=[x0(1) x(1)];
x2=[x0(2) x(2)];
x3=[x0(3) x(3)];
y1=[y0(1) y(1)];
y2=[y0(2) y(2)];
y3=[y0(3) y(3)];
line(x1,y1,...
    'LineWidth',6,...
    'Color',color);
hold on;
line(x2,y2,...
    'LineWidth',6,...
    'Color',color);
hold on;
line(x3,y3,...
    'LineWidth',6,...
    'Color',color);
hold on;
title('relative error on fluct', 'fontsize', 12);
xlabel('relative error', 'fontsize', 12);
ylabel('%', 'fontsize', 12);
xlim([0,4])
hold off;
saveas(d, [folders 'fluct.png']);
close

m=figure;
x = M(:,1);
y = M(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(m, [folders 'error(eps).png']);
close

n=figure;
x = N(:,1);
y = N(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('iter on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('iter', 'fontsize', 12);
hold off;
saveas(n, [folders 'iter(eps).png']);
close

o=figure;
x = O(:,1);
y = O(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('error(iter)', 'fontsize', 12);
xlabel('iter', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(o, [folders 'error(iter).png']);
close

oo=figure;
x = OO(:,1);
y = OO(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',color,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color);
title('error(p)', 'fontsize', 12);
xlabel('p', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(oo, [folders 'error(p).png']);
close