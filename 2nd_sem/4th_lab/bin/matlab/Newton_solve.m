function x = Newton_solve(F, J, x0)
    max_iterations = 100;
    tolerance = 1e-6;
    x = x0;
    
    for k = 1:max_iterations
        delta_x = -J \ F(x);
        x = x + delta_x;
        
        if norm(delta_x) < tolerance
            break;
        end
    end
end
