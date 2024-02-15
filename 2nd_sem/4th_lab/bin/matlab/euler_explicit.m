function [Y, T] = euler_explicit(f, t0, b, y0, h)
    T = t0 : h : b;
    Y = zeros(size(T));
    Y(1) = y0;
    for n = 1 : length(T) - 1
        Y(n+1) = Y(n) + h * f(T(n), Y(n));
    end
end
