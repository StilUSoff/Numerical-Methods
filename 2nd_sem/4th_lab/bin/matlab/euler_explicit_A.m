function result = euler_explicit_A(A, T, x0)
    result = zeros(size(T, 2), size(x0, 1));
    result(1, :) = x0;
    for i = 2:size(T, 2)
        result(i, :) = result(i - 1, :)' + (T(i)-T(i-1)) * (A * result(i - 1, :)');
    end
end
