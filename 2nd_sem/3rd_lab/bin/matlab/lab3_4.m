color1=['#bbafff'];
color2=['#afffe3'];
color3=['#c9ffaf'];

fun = @(x) x.^5-5.2.*x.^3+2.5.*x.^2-7.*x
fun_2 = @(x,y) y.^3-5.2.*x.^3+2.5.*x.^2-7.*x
fun_5 = @(x,y,z) y.^3-5.2.*x.^3+2.5.*x.^2-7.*z.^2
fun_6 = @(x) 1./(x.^2+x-2.0)
fun_sin = @(x) sin(x.^2).*(x.^5-5.2.*x.^3+2.5.*x.^2-7.*x)
fun_abs = @(x) abs(x.^5-5.2.*x.^3+2.5.*x.^2-7.*x)





a=-2.1;
b=2.1;
real=3087/200;
real_2=52.84635300000001;
real_3=2.596124661372699;
real_31=2.596124661372699;
real_4=64827/1000;
real_5=-12252303/25000;
real_6=log(4)/3.0

%1

for p=2:10
    integ=integral(fun,a,b,'AbsTol',10^(-p));
    error(p-1)=(abs(real-integ));
    eps(p-1)=10^(-p);
end

for p=2:10
    abs_integ=integral(fun_abs,a,b,'AbsTol',10^(-p), 'RelTol',10^(-p));
    abs_error(p-1)=(abs(real_2-abs_integ));
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
hold on;
loglog(eps,abs_error,'-*',...
    'LineWidth',2,...
    'Color',color2,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color2);
hold on;
legend('P(x)', '|P(x)|');
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
hold off;
saveas(f, ['results/1st.png']);
close



%2

for p=2:10
    integ=integral(fun_sin,a,b,'AbsTol',10^(-p),'RelTol', 10^(-p));
    error(p-1)=(abs(real_3-integ));
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
saveas(f, ['results/2nd.png']);
close


%3

for p=2:10
    integ=quadgk(fun_sin,a,b,'RelTol',10^(-p),'AbsTol',10^(-p));
    error(p-1)=(abs(real_31-integ));
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
saveas(f, ['results/3rd.png']);
close


%4

for p=2:10
    integ=integral2(fun_2,a,b,a,b,'RelTol',10^(-p),'AbsTol',10^(-p));
    error(p-1)=(abs(real_4-integ));
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
saveas(f, ['results/4th.png']);
close


%5

for p=2:10
    integ=integral3(fun_5,a,b,a,b,a,b,'RelTol',10^(-p),'AbsTol',10^(-p));
    error(p-1)=(abs(real_5-integ));
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
saveas(f, ['results/5th.png']);
close


%6

a = 2;
b = 200;
for p=2:10
    integ=integral(fun_6,a,b,'AbsTol',10^(-p));
    error(p-1)=(abs(real_6-integ));
    eps(p-1)=10^(-p);
end

f=figure;
loglog(eps,error,'-*',...
    'LineWidth',2,...
    'Color',color1,...
    'MarkerSize',4,... 
    'MarkerFaceColor',color1);
title('error on eps', 'fontsize', 12);
xlabel('eps', 'fontsize', 12);
ylabel('error', 'fontsize', 12);
saveas(f, ['results/6th.png']);
close