function [Y, T] = euler_implicit(f, t0, b, y0, h)
    T = t0 : h : b;
    Y = zeros(size(T));
    Y(1) = y0;
    for n = 1 : length(T) - 1
        f1 = @(x) x - Y(n) - h * f(T(n+1), x);
        Y(n+1) = fsolve(f1, 0);
    end
end

