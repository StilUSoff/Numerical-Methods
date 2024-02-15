!folder = ['bin/researches_res/'];
folder = ['researches_res/'];
!folders = ['bin/graphics/'];
folders = ['graphics/'];
        
J = load([folder 'error.txt']);
K = load([folder 'time.txt']);
L = load([folder 'Afluct.txt']);   
D = load([folder 'Bfluct.txt']);

j=figure;
x = J(:,1);
y = J(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',"#c5afff",...
    'MarkerSize',4,... 
    'MarkerFaceColor',"#c5afff");
title('error(condition number)', 'fontsize', 12);
xlabel('condition number', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(j, [folders 'error.png']);
close
    
k=figure;
x = K(:,1);
y = K(:,2);
loglog(x,y,'-*',...
    'LineWidth',4,...
    'Color',"#c5afff",...
    'MarkerSize',4,... 
    'MarkerFaceColor',"#c5afff");
title('time(condition number)', 'fontsize', 12);
xlabel('condition number', 'fontsize', 12);
ylabel('time', 'fontsize', 12);
hold off;
saveas(k, [folders 'time.png']);
close

l=figure;
x = L(:,1)';
y = L(:,2)';
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
    'Color',"#c5afff");
hold on;
line(x2,y2,...
    'LineWidth',6,...
    'Color',"#c5afff");
hold on;
line(x3,y3,...
    'LineWidth',6,...
    'Color',"#c5afff");
hold on;
title('fluct from A', 'fontsize', 12);
xlabel('fluct', 'fontsize', 12);
ylabel('%', 'fontsize', 12);
xlim([0,4])
hold off;
saveas(l, [folders 'Afluct.png']);
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
    'Color',"#c5afff");
hold on;
line(x2,y2,...
    'LineWidth',6,...
    'Color',"#c5afff");
hold on;
line(x3,y3,...
    'LineWidth',6,...
    'Color',"#c5afff");
hold on;
title('fluct from B', 'fontsize', 12);
xlabel('fluct', 'fontsize', 12);
ylabel('%', 'fontsize', 12);
xlim([0,4])
hold off;
saveas(d, [folders 'Bfluct.png']);
close





