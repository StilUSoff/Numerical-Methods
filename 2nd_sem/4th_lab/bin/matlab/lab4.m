function lab4_matlab()
format longg
for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode45(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
end
f1 = plot_graphics(log10(Eps), log10(error), 'ode45', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode113(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
end
f2 = plot_graphics(log10(Eps), log10(error), 'ode113', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode15s(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
%     f(i) = figure;
%     plot(x,y, '-o')
%     hold on
%     plot(x1, y1, '-o')
%     hold off
end
f4 = plot_graphics(log10(Eps), log10(error), 'ode15s', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode23(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
%     f(i) = figure;
%     plot(x,y, '-o')
%     hold on
%     plot(x1, y1, '-o')
%     hold off
end
f5 = plot_graphics(log10(Eps), log10(error), 'ode23', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode23s(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
%     f(i) = figure;
%     plot(x,y, '-o')
%     hold on
%     plot(x1, y1, '-o')
%     hold off
end
f6 = plot_graphics(log10(Eps), log10(error), 'ode23s', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode23t(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
%     f(i) = figure;
%     plot(x,y, '-o')
%     hold on
%     plot(x1, y1, '-o')
%     hold off
end
f7 = plot_graphics(log10(Eps), log10(error), 'ode23t', 'log10(Epsilon)', 'log10(Error)')

for i = 1:16
    Eps(i) = 10^-i;
    options = odeset('RelTol', Eps(i));
    tspan = [1,2];
    y0 = 1;
    [x, y] = ode23tb(@(x,y) odu(x,y), tspan, y0, options);
    y1 = true_solution(x);
    error(i) = max(abs(y-y1));
%     f(i) = figure;
%     plot(x,y, '-o')
%     hold on
%     plot(x1, y1, '-o')
%     hold off
end
f8 = plot_graphics(log10(Eps), log10(error), 'ode23tb', 'log10(Epsilon)', 'log10(Error)')


j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
tt = 0:0.005:2;
yy = true_solution(tt);
plot(tt, yy)
hold on

a = 0;
b = 2;
h = 0.001;
y0 = 1

f = @(t, y) odu(t, y);

[Y, T] = euler_explicit(f, a, b, y0, h);
plot(T, Y, 'r')
title('Explicit Eueler')
hold off
saveas(j, ['graphics/Explicit Eueler.png'])
close


f = @(t, y) 5*(y - t^2);
yt = @(t) t.^2 + 0.4*t + 0.08;
t0 = 1;
b = 2;
y0 = 0.08;

tt = t0:0.005:b;
yy = yt(tt);
j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(tt, yy)
title('Explicit Eueler step = 0.08:-0.01:0.01')
legend('true', 'Location', 'best')
hold on

i=1;
for h = 0.08:-0.01:0.01
    [Y, T] = euler_explicit(f, t0, b, y0, h);
	yy = yt(T);
    h_step1(i) = h;
    error_1(i) = max(abs(Y-yy));
    plot(T, Y, 'DisplayName', ['h=',num2str(h)])
    i=i+1;
end
saveas(j, ['graphics/Explicit Eueler step = 0.08:-0.01:0.01.png'])
close


f = @(t, y) 5*(y - t^2);
yt = @(t) t.^2 + 0.4*t + 0.08;
t0 = 0.08;
b = 0.5;
y0 = 0.08;

tt = t0:0.005:b;
yy = yt(tt);
j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(tt, yy)
title('Explicit Eueler different steps on 0.8:5')
legend('true', 'Location', 'best')
hold on

h = 0.1;
for i = 1:8
    h = h / 2;
    [Y, T] = euler_explicit(f, t0, b, y0, h);
    yy = yt(T);
    error_2(i) = max(abs(Y-yy));
    h_step2(i) = h;
    plot(T, Y, 'DisplayName', ['h=',num2str(h)])
end
saveas(j, ['graphics/Explicit Eueler different steps on 0.8:5.png'])
close


j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(h_step1, log10(error_1))
title('error on step')
saveas(j, ['graphics/Explicit error on step.png'])
close

j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(h_step2, log10(error_2))
title('error on step 2')
saveas(j, ['graphics/Explicit error on step 2.png'])
close


j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
tt = 0:0.005:2;
yy = true_solution(tt);
plot(tt, yy)
hold on

a = 0;
b = 2;
h = 0.001;
y0 = 1;

f = @(t, y) odu(t, y);

[Y, T] = euler_implicit(f, a, b, y0, h);
plot(T, Y, 'r')
title('Implicit Eueler')
hold off
saveas(j, ['graphics/Implicit Eueler.png'])
close


f = @(t, y) -100 * y + 10;
yt = @(t) 0.1 + 0.9 * exp(-100 * t);
t0 = 0;
b = 10;
y0 = 1;

tt = t0:0.005:b;
yy = yt(tt);
j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(tt, yy)
title('Implicit Eueler step = 0.08:-0.01:0.01')
legend('true', 'Location', 'best')
hold on

i=1;
for h = 0.08:-0.01:0.01
    [Y, T] = euler_explicit(f, t0, b, y0, h);
	yy = yt(T);
    h_step1(i) = h;
    error_1(i) = max(abs(Y-yy));
    plot(T, Y, 'DisplayName', ['h=',num2str(h)])
    i=i+1;
end
saveas(j, ['graphics/Implicit Eueler step = 0.08:-0.01:0.01.png'])
close

f = @(t, y) 5*(y - t^2);
yt = @(t) t.^2 + 0.4*t + 0.08;
t0 = 0.08;
b = 1;
y0 = 1;

tt = t0:0.005:b;
yy = yt(tt);
j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(tt, yy)
title('Implicit Eueler different steps on 0.8:5')
legend('true', 'Location', 'best')
hold on

h = 0.1;
for i = 1:8
    h = h / 2;
    [Y, T] = euler_explicit(f, t0, b, y0, h);
    yy = yt(T);
    error_2(i) = max(abs(Y-yy));
    h_step2(i) = h;
    plot(T, Y, 'DisplayName', ['h=',num2str(h)])
end
saveas(j, ['graphics/Implicit Eueler different steps on 0.8:5.png'])
close

j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(h_step1, log10(error_1))
title('error on step')
saveas(j, ['graphics/Implicit error on step.png'])
close

j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(h_step2, log10(error_2))
title('error on step 2')
saveas(j, ['graphics/Implicit error on step 2.png'])
close

k = 1;
A = rand(10, 10);
for i = 2 : 0.2105 : 6
    lamda = linspace(-1, round(-10^i), 10);
    max_length(k) = 10.0^-i;
    [Q, R] = qr(A);
    R = R - diag(diag(R)) + diag(lamda);
    A = Q*R*Q';
    x0 = randi([0,1],10,1);
    a = 1;
    b = 2;
    [t, true_sol_ode45] = ode45(@(t, x) A_odu(t, x, A), [a,b] , x0);
    solution = euler_explicit_A(A, t', x0);
    solution2 = euler_implicit_A(A, t', x0);
    error1(k) = max(norm(abs(solution - true_sol_ode45)));
    error2(k) = max(norm(abs(solution2 - true_sol_ode45)));
    k = k + 1;
end
j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(max_length, log10(error1))
title('error explicit A')
saveas(j, ['graphics/Explicit error A.png'])
close


j=figure('position', [350 350 1550 950],'DefaultAxesFontSize',14);
plot(max_length, log10(error2))
title('error implicit A')
saveas(j, ['graphics/Implicit error A.png'])
close

end