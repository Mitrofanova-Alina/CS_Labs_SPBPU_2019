function [R] = find_curv_radius(A, B, C)
    a = B - A;
    b = C - B;
    c = C - A;
    R = norm(a) / (2 * sqrt(1 - (sum((c) .* (b)) / (norm(c) * norm(b)))^2));
end