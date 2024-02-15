function result = euler_implicit_A(A, T, x0)
    result = zeros(size(T, 2), size(x0, 1));
    x = x0;
    result(1, :) = x;
    for i = 2:length(T)
        F = @(x_next) x_next - x - (T(i)-T(i-1)) * (A * x_next);
        J = eye(length(x)) - (T(i)-T(i-1)) * A;
        x_next = Newton_solve(F, J, x);
        
        x = x_next;
        result(i, :) = x;
    end
end

