folder = ['researches_res/'];
folders = ['graphics/'];
color1=['#B3FFE6'];
color2=['#ffdfb3'];
color3=['#E6B3FF'];
color4=['#bffa91'];
color5=['#FE3A3F'];
color6=['#29ACB4'];


Aepsabs = load([folder '1_err_on_eps_abs.txt']);
Aeps = load([folder '1_err_on_eps.txt']);
Alen = load([folder '1_err_on_len.txt']);
Aiterabs = load([folder '1_iter_on_eps_abs.txt']);
Aiter = load([folder '1_iter_on_eps.txt']);
B1er = load([folder '2_1_err_on_eps.txt']);
B1iter = load([folder '2_1_iter_on_eps.txt']);
B1len = load([folder '2_1_len_on_x.txt']);
B2erfl1 = load([folder '2_2_err_fluct_1.txt']);
B2erfl2 = load([folder '2_2_err_fluct_2.txt']);
B2erfl3 = load([folder '2_2_err_fluct_3.txt']);
B2er = load([folder '2_2_err_on_eps.txt']);
B2iter = load([folder '2_2_iter_on_eps.txt']);
B2lenfl1 = load([folder '2_2_len_fluct_1.txt']);
B2lenfl2 = load([folder '2_2_len_fluct_2.txt']);
B2lenfl3 = load([folder '2_2_len_fluct_3.txt']);
B2len = load([folder '2_2_len_on_x.txt']);
B3er = load([folder '2_3_err_on_eps.txt']);
B3iter = load([folder '2_3_iter_on_eps.txt']);
B3len = load([folder '2_3_len_on_x.txt']);
C1ereps = load([folder '3_3_err_on_eps_1.txt']);
C1erepsabs = load([folder '3_3_err_on_eps_2.txt']);
C1richard = load([folder '3_3_err_on_eps_richard.txt']);
C1iter = load([folder '3_3_iter_on_eps__1.txt']);
C1iterabs = load([folder '3_3_iter_on_eps__2.txt']);
C2ereps = load([folder '3_4_err_on_eps_1.txt']);
C2erepsabs = load([folder '3_4_err_on_eps_2.txt']);
C2richard = load([folder '3_4_err_on_eps_richard.txt']);
C2iter = load([folder '3_4_iter_on_eps__1.txt']);
C2iterabs = load([folder '3_4_iter_on_eps__2.txt']);


a1=figure;
x = Aeps(:,2);
y = Aeps(:,1);
xx = Aepsabs(:,2);
yy = Aepsabs(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
legend('error f(x)', 'error |f(x)|','Location','northwest');
title('error on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(error)', 'fontsize', 12);
hold off;
saveas(a1, [folders '1_error_eps.png']);
close
 

a2=figure;
x = Aiter(:,2);
y = Aiter(:,1);
xx = Aiterabs(:,2);
yy = Aiterabs(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
legend('iter f(x)', 'iter |f(x)|');
title('iter on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(iter)', 'fontsize', 12);
hold off;
saveas(a2, [folders '1_iter_eps.png']);
close


a3=figure;
x = Alen(:,2);
y = Alen(:,1);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
title('error on lenght', 'fontsize', 12);
xlabel('lenght', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(a3, [folders '1_error_len.png']);
close


b1=figure;
x = B1er(:,2);
y = B1er(:,1);
xx = B2er(:,2);
yy = B2er(:,1);
xxx = B3er(:,2);
yyy = B3er(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
plot(log2(xxx),log2(yyy),'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
legend('error f(x)', 'error I integral','error II integral','Location','northwest');
title('error on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(error)', 'fontsize', 12);
hold off;
saveas(b1, [folders '2_error_eps.png']);
close


b2=figure;
x = B1iter(:,2);
y = B1iter(:,1);
xx = B2iter(:,2);
yy = B2iter(:,1);
xxx = B3iter(:,2);
yyy = B3iter(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
plot(log2(xxx),log2(yyy),'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
legend('iter f(x)', 'iter I integral','iter II integral');
title('iter on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(iter)', 'fontsize', 12);
hold off;
saveas(b2, [folders '2_iter_eps.png']);
close


b3_1=figure;
x = B1len(:,1);
y = B1len(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
title('length f(x)', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b3_1, [folders '2_1_len_x.png']);
close

b3_2=figure;
x = B2len(:,1);
y = B2len(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
title('length I integral', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b3_2, [folders '2_2_len_x.png']);
close

b3_3=figure;
x = B3len(:,1);
y = B3len(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
title('length II integral', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b3_3, [folders '2_3_len_x.png']);
close


b4=figure;
x = B2erfl1(:,2);
y = B2erfl1(:,1);
xx = B2erfl2(:,2);
yy = B2erfl2(:,1);
xxx = B2erfl3(:,2);
yyy = B2erfl3(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
plot(log2(xxx),log2(yyy),'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
legend('error 1%', 'error 2%','error 3%','Location','northwest');
title('error on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(error)', 'fontsize', 12);
hold off;
saveas(b4, [folders '2_fluct_error_eps.png']);
close


b5_1=figure;
x = B2lenfl1(:,1);
y = B2lenfl1(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
title('length 1%', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b5_1, [folders '2_fluct_1_len_x.png']);
close

b5_2=figure;
x = B2lenfl2(:,1);
y = B2lenfl2(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
title('length 2%', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b5_2, [folders '2_fluct_2_len_x.png']);
close

b5_3=figure;
x = B2lenfl3(:,1);
y = B2lenfl3(:,2);
loglog(x,y,'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
title('length 3%', 'fontsize', 12);
xlabel('x', 'fontsize', 12);
ylabel('length', 'fontsize', 12);
hold off;
saveas(b5_3, [folders '2_fluct_3_len_x.png']);
close


c1=figure;
x = C1ereps(:,2);
y = C1ereps(:,1);
xx = C2ereps(:,2);
yy = C2ereps(:,1);
xxx = C1erepsabs(:,2);
yyy = C1erepsabs(:,1);
xxxx = C2erepsabs(:,2);
yyyy = C2erepsabs(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
plot(log2(xxx),log2(yyy),'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
plot(log2(xxxx),log2(yyyy),'-*',...
    'LineWidth',2,...
    'Color',color4,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color4);
hold on;
legend('error n=3', 'error n=4','error n=3 (abs)','error n=4 (abs)','Location','northwest');
title('error on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(error)', 'fontsize', 12);
hold off;
saveas(c1, [folders '3_error_eps.png']);
close


c2=figure;
x = C1richard(:,2);
y = C1richard(:,1);
xx = C2richard(:,2);
yy = C2richard(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
legend('error n=3', 'error n=4','Location','northwest');
title('Richardson error on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(error)', 'fontsize', 12);
hold off;
saveas(c2, [folders '3_richardson.png']);
close


c3=figure;
x = C1iter(:,2);
y = C1iter(:,1);
xx = C2iter(:,2);
yy = C2iter(:,1);
xxx = C1iterabs(:,2);
yyy = C1iterabs(:,1);
xxxx = C2iterabs(:,2);
yyyy = C2iterabs(:,1);
plot(log2(x),log2(y),'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
plot(log2(xx),log2(yy),'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
plot(log2(xxx),log2(yyy),'-*',...
    'LineWidth',2,...
    'Color',color3,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color3);
hold on;
plot(log2(xxxx),log2(yyyy),'-*',...
    'LineWidth',2,...
    'Color',color4,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color4);
hold on;
legend('iter n=3', 'iter n=4','iter n=3 (abs)','iter n=4 (abs)','Location','northeast');
title('iter on eps', 'fontsize', 12);
xlabel('log_2(eps)', 'fontsize', 12);
ylabel('log_2(iter)', 'fontsize', 12);
hold off;
saveas(c3, [folders '3_iter_eps.png']);
close