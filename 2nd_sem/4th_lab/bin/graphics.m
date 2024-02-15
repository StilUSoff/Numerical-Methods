%folder = ['bin/researches_res/'];
folder = ['researches_res/'];
%folders = ['bin/graphics/'];
folders = ['graphics/'];
color=['#e6b3ff'];
color2=['#B3FFE6'];
color3=['#000000'];

% J = load([folder '1_t_y_0_1.txt']);
% K = load([folder '1_t_y_0_05.txt']);
% D = load([folder '2_tcoord_step.txt']);
M = load([folder '3_eps_error.txt']);

x_r = 1 : 2.4/100000.0 : 2.0 ;
y_r = 1.74*exp(x_r.^2)-x_r.^2-1;

j=figure;
x = J(:,1);
y = J(:,2);
x2 = K(:,1);
y2 = K(:,2);
loglog(x,y,'-',...
    'LineWidth',5,...
    'Color',color);
hold on;
loglog(x2,y2,'-',...
    'LineWidth',3,...
    'Color',color2);
hold on;
loglog(x_r,y_r,'-',...
    'LineWidth',1,...
    'Color',color3);
title('y(x)', 'fontsize', 12);
legend('h=10^-^2','h=10^-^6','true solution');
xlabel('x', 'fontsize', 12);
ylabel('y', 'fontsize', 12);
hold off;
saveas(j, [folders '1_t_y.png']);
close
    
d=figure;
x = D(:,1);
y = D(:,2);
loglog(x,y,'-',...
    'LineWidth',1,...
    'Color',color);
title('lenght(x)', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('lenght', 'fontsize', 12);
hold off;
saveas(d, [folders '2_tcoord_step.png']);
close

m=figure;
x = M(:,1);
y = abs(M(:,2)-1.74*exp(2.^2)+2.^2+1);
loglog(x,y,'-|',...
    'LineWidth',1,...
    'Color',color);
title('error(eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error, c', 'fontsize', 12);
hold off;
saveas(m, [folders '3_eps_error.png']);
close